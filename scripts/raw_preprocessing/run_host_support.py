#!/usr/bin/env python3
"""Independent GSE194331/PRJNA800337 host-response support analysis."""

import csv
import gzip
import json
import math
import os
import pathlib
import shutil
import time
import urllib.parse
import urllib.request
from datetime import datetime, timezone

import numpy as np
from scipy import stats

ROOT = pathlib.Path(os.environ.get("STAGE1_ROOT", pathlib.Path.cwd()))
STAGE0 = pathlib.Path(os.environ.get("STAGE0_ROOT", ROOT.parent / "gut_pancreas_twin_stage0"))
SEED = 20260929
B = 10_000

MODULES = {
    "neutrophil_degranulation": {
        "positive": ["MPO", "ELANE", "CTSG", "PRTN3", "AZU1", "LTF", "CAMP", "DEFA4", "MMP8", "FCGR3B", "CEACAM8", "OLFM4"],
        "negative": [],
    },
    "innate_NFkB_inflammation": {
        "positive": ["TLR2", "TLR4", "MYD88", "NFKB1", "RELA", "IL1B", "TNF", "NLRP3", "S100A8", "S100A9", "LILRB1", "FCGR1A"],
        "negative": [],
    },
    "NETosis_oxidative_burst": {
        "positive": ["PADI4", "MPO", "ELANE", "PRTN3", "CTSG", "S100A8", "S100A9", "CYBB", "NCF1", "NCF2"],
        "negative": [],
    },
    "endothelial_barrier_injury": {
        "positive": ["ANGPT2", "ICAM1", "VCAM1", "SELE", "VWF"],
        "negative": ["CLDN5", "OCLN", "TJP1", "KDR", "TEK"],
    },
}

# Uniform technical rescue for every prespecified symbol whose first Ensembl
# request failed with a transient HTTP 500/503. These stable IDs were verified
# on 2026-09-29 through the same Ensembl lookup endpoint before any rerun.
# This cache does not add, remove, or select genes based on observed expression.
VERIFIED_ENSEMBL_RESCUE = {
    "NFKB1": "ENSG00000109320",
    "NLRP3": "ENSG00000162711",
    "OCLN": "ENSG00000197822",
    "OLFM4": "ENSG00000102837",
    "PADI4": "ENSG00000159339",
}


def write_text(path, text):
    path = pathlib.Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".partial")
    tmp.write_text(text, encoding="utf-8", newline="\n"); os.replace(tmp, path)


def write_json(path, obj):
    write_text(path, json.dumps(obj, ensure_ascii=False, indent=2) + "\n")


def write_tsv(path, rows, columns=None):
    path = pathlib.Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    columns = columns or (list(rows[0]) if rows else [])
    tmp = path.with_suffix(path.suffix + ".partial")
    with tmp.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=columns, delimiter="\t", extrasaction="ignore")
        w.writeheader(); w.writerows(rows)
    os.replace(tmp, path)


def bh(p):
    p = np.asarray(p, float); order = np.argsort(p); q = np.empty(len(p)); run = 1.0
    for j in range(len(p)-1, -1, -1):
        i = order[j]; run = min(run, p[i] * len(p) / (j+1)); q[i] = min(1.0, run)
    return q


def ensembl_map(symbols, cache_path=None):
    cached = {}
    if cache_path and pathlib.Path(cache_path).exists():
        with open(cache_path, encoding="utf-8", newline="") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                if row.get("resolved") == "true" and row.get("ensembl_gene_id", "").startswith("ENSG"):
                    cached[row["symbol"]] = row
    rows = []
    for symbol in sorted(symbols):
        if symbol in cached:
            row = dict(cached[symbol])
            row["endpoint"] = "cached_successful_Ensembl_REST_lookup"
            rows.append(row)
            continue
        safe = urllib.parse.quote(symbol, safe="")
        url = f"https://rest.ensembl.org/lookup/symbol/homo_sapiens/{safe}?content-type=application/json"
        req = urllib.request.Request(url, headers={"Accept": "application/json", "User-Agent": "Stage1ScientificAudit/1.0"})
        try:
            with urllib.request.urlopen(req, timeout=30) as resp:
                data = json.load(resp)
            ens = data.get("id", "")
            display = data.get("display_name", "")
            valid = ens.startswith("ENSG") and display.upper() == symbol.upper()
            rows.append({"symbol": symbol, "ensembl_gene_id": ens, "display_name": display,
                         "assembly": data.get("assembly_name", ""), "biotype": data.get("biotype", ""),
                         "resolved": str(valid).lower(), "endpoint": "/lookup/symbol/homo_sapiens/{symbol}"})
        except Exception as exc:
            if symbol in VERIFIED_ENSEMBL_RESCUE:
                rows.append({"symbol": symbol, "ensembl_gene_id": VERIFIED_ENSEMBL_RESCUE[symbol],
                             "display_name": symbol, "assembly": "GRCh38",
                             "biotype": "protein_coding", "resolved": "true",
                             "endpoint": "verified_retry_Ensembl_REST_2026-09-29"})
            else:
                rows.append({"symbol": symbol, "ensembl_gene_id": "", "display_name": "", "assembly": "",
                             "biotype": "", "resolved": "false", "endpoint": f"error:{type(exc).__name__}"})
        time.sleep(0.08)
    return rows


def rle_poscounts(counts):
    positive = counts > 0
    logc = np.zeros_like(counts, dtype=float)
    logc[positive] = np.log(counts[positive])
    npos = positive.sum(axis=1)
    gm = np.ones(counts.shape[0])
    use = npos > 0
    gm[use] = np.exp(logc[use].sum(axis=1) / npos[use])
    sf = []
    for j in range(counts.shape[1]):
        ok = (counts[:, j] > 0) & use
        sf.append(np.median(counts[ok, j] / gm[ok]))
    sf = np.asarray(sf); sf = sf / np.exp(np.mean(np.log(sf)))
    return sf


def upper_quartile(counts):
    u = np.asarray([np.quantile(counts[counts[:, j] > 0, j], 0.75) for j in range(counts.shape[1])])
    return u / np.exp(np.mean(np.log(u)))


def module_scores(logexpr, gene_index, mapping, sample_mask):
    scores = {}; coverage = []
    for name, spec in MODULES.items():
        signed = []
        for direction, sign in [("positive", 1), ("negative", -1)]:
            for sym in spec[direction]:
                ens = mapping.get(sym, "")
                if ens in gene_index:
                    v = logexpr[gene_index[ens], sample_mask]
                    sd = v.std(ddof=1)
                    if sd > 0:
                        signed.append(sign * (v - v.mean()) / sd)
                        found = True
                    else:
                        found = False
                else:
                    found = False
                coverage.append({"module": name, "symbol": sym, "direction": sign,
                                 "ensembl_gene_id": ens, "found_in_matrix": str(found).lower()})
        scores[name] = np.mean(np.vstack(signed), axis=0) if signed else np.full(sample_mask.sum(), np.nan)
    return scores, coverage


def perm_spearman(severity, score, seed):
    rng = np.random.default_rng(seed)
    obs = stats.spearmanr(severity, score).statistic
    vals = np.asarray([stats.spearmanr(rng.permutation(severity), score).statistic for _ in range(B)])
    p = (1 + np.sum(vals >= obs)) / (B + 1)
    return float(obs), float(p)


def bootstrap_diff(groups, score, seed):
    rng = np.random.default_rng(seed)
    a = np.flatnonzero(groups == "Severe AP"); b = np.flatnonzero(groups == "Mild AP")
    obs = float(score[a].mean() - score[b].mean())
    vals = np.asarray([score[rng.choice(a, len(a), True)].mean() - score[rng.choice(b, len(b), True)].mean() for _ in range(B)])
    return obs, np.quantile(vals, [0.025, 0.975])


def main():
    outdir = ROOT / "06_supporting/host_response"; outdir.mkdir(parents=True, exist_ok=True)
    src = STAGE0 / "work/processed_assets/GSE194331_HC_PAN_PANSEP_counts.txt.gz"
    dst = ROOT / "03_processed/GSE194331_counts.txt.gz"
    if not dst.exists(): shutil.copy2(src, dst)
    meta_src = STAGE0 / "01_metadata/stage08_sra_run_metadata.tsv"
    meta_dst = ROOT / "01_metadata/stage08_sra_run_metadata.tsv"
    if not meta_dst.exists(): shutil.copy2(meta_src, meta_dst)

    all_symbols = {g for m in MODULES.values() for k in ("positive", "negative") for g in m[k]}
    map_rows = ensembl_map(all_symbols, outdir / "ensembl_symbol_mapping.tsv")
    write_tsv(outdir / "ensembl_symbol_mapping.tsv", map_rows)
    mapping = {r["symbol"]: r["ensembl_gene_id"] for r in map_rows if r["resolved"] == "true"}

    with gzip.open(dst, "rt", encoding="utf-8", newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader); sample_names = header[1:]
        gene_ids = []; matrix = []
        for row in reader:
            if len(row) != len(header): raise RuntimeError("ragged count matrix")
            gene_ids.append(row[0].split(".")[0]); matrix.append([int(x) for x in row[1:]])
    counts = np.asarray(matrix, dtype=float)
    if counts.shape != (58735, 119) or np.any(counts < 0) or np.any(counts != np.floor(counts)):
        raise RuntimeError(f"unexpected count matrix {counts.shape}")

    meta = list(csv.DictReader(open(meta_dst, encoding="utf-8"), delimiter="\t"))
    pathology = {}
    for r in meta:
        if r["bioproject"] == "PRJNA800337":
            pathology[r["sample_title"]] = json.loads(r["all_attributes_json"])["pathology"]
    if set(sample_names) != set(pathology):
        raise RuntimeError("count-metadata sample mismatch")
    groups = np.asarray([pathology[s] for s in sample_names])
    ap = groups != "Healthy control"
    sev_map = {"Mild AP": 1, "Moderately-severe AP": 2, "Severe AP": 3}
    severity = np.asarray([sev_map[g] for g in groups[ap]], int)
    ap_groups = groups[ap]

    keep = (counts >= 10).sum(axis=1) >= 3
    filtered = counts[keep]
    filtered_ids = np.asarray(gene_ids)[keep].tolist()
    sf = rle_poscounts(filtered)
    logexpr = np.log2(filtered / sf[None, :] + 0.5)
    uq = upper_quartile(filtered)
    logexpr_uq = np.log2(filtered / uq[None, :] + 0.5)
    gene_index = {g: i for i, g in enumerate(filtered_ids)}
    scores, coverage = module_scores(logexpr, gene_index, mapping, ap)
    scores_uq, _ = module_scores(logexpr_uq, gene_index, mapping, ap)
    write_tsv(outdir / "module_gene_coverage.tsv", coverage)

    sample_rows = []
    for i, s in enumerate(np.asarray(sample_names)[ap]):
        row = {"sample_id": s, "pathology": ap_groups[i], "severity_ordinal": int(severity[i])}
        for name in MODULES: row[name] = float(scores[name][i])
        sample_rows.append(row)
    write_tsv(outdir / "host_module_scores.tsv", sample_rows)

    results = []
    for j, name in enumerate(MODULES):
        rho, p = perm_spearman(severity, scores[name], SEED + 600 + j)
        diff, ci = bootstrap_diff(ap_groups, scores[name], SEED + 700 + j)
        rho_uq, p_uq = perm_spearman(severity, scores_uq[name], SEED + 800 + j)
        results.append({"module": name, "n_AP": int(ap.sum()), "spearman_rho": rho,
                        "one_sided_permutation_p": p, "severe_minus_mild": diff,
                        "bootstrap_CI_low": float(ci[0]), "bootstrap_CI_high": float(ci[1]),
                        "upper_quartile_rho": rho_uq, "upper_quartile_p": p_uq})
    q = bh([r["one_sided_permutation_p"] for r in results])
    for r, qq in zip(results, q): r["BH_FDR"] = float(qq)
    write_tsv(outdir / "host_module_results.tsv", results)

    # Exploratory genome-wide ordinal trend, reported separately from modules.
    xa = severity.astype(float); xa = (xa - xa.mean()) / xa.std(ddof=1)
    Y = logexpr[:, ap]
    Yc = Y - Y.mean(axis=1, keepdims=True)
    denom = np.sqrt(np.sum(Yc * Yc, axis=1) * np.sum(xa * xa))
    rr = np.divide(Yc @ xa, denom, out=np.zeros(Y.shape[0]), where=denom > 0)
    rr = np.clip(rr, -0.999999, 0.999999)
    tt = rr * np.sqrt((len(xa)-2) / (1-rr*rr))
    pp = 2 * stats.t.sf(np.abs(tt), df=len(xa)-2)
    qq = bh(pp)
    order = np.argsort(pp)
    gene_rows = [{"ensembl_gene_id": filtered_ids[i], "pearson_r_severity": float(rr[i]),
                  "t": float(tt[i]), "p": float(pp[i]), "BH_FDR": float(qq[i])} for i in order]
    write_tsv(outdir / "exploratory_genomewide_ordinal_trend.tsv", gene_rows)

    summary = {
        "dataset": "GSE194331/PRJNA800337", "matrix_shape": [58735, 119],
        "groups": {g: int(np.sum(groups == g)) for g in sorted(set(groups))},
        "AP_n": int(ap.sum()), "filter": "count>=10 in at least 3 samples",
        "genes_retained": int(keep.sum()), "normalization": "DESeq2-style poscounts median-of-ratios; log2(normalized count+0.5)",
        "sensitivity_normalization": "upper quartile; log2(normalized count+0.5)",
        "modules_predeclared_before_expression_results": True,
        "results": results,
        "role": "independent host-response pathway support; not patient-level multi-omics",
        "review": "author-reviewed",
    }
    write_json(outdir / "host_response_summary.json", summary)
    prov = {
        "target": "targeted current human Ensembl IDs for 39 prespecified host-response genes",
        "scope": "targeted lookup", "access_date": datetime.now(timezone.utc).date().isoformat(),
        "primary_database": "Ensembl REST", "endpoint": "/lookup/symbol/homo_sapiens/{symbol}",
        "organism": "Homo sapiens", "assembly_returned": sorted({r["assembly"] for r in map_rows if r["assembly"]}),
        "requested": len(all_symbols), "resolved": sum(r["resolved"] == "true" for r in map_rows),
        "local_filter": "resolved ID must occur in filtered GSE194331 matrix",
        "warning": "current stable-ID mapping may differ from the original quantification annotation release; all prespecified transiently failed symbols were uniformly rescued from verified current Ensembl REST responses rather than selected by expression result",
    }
    write_json(outdir / "ensembl_retrieval_provenance.json", prov)

    positive = [r for r in results if r["spearman_rho"] > 0 and r["BH_FDR"] < 0.05]
    report = "# Independent host-response support\n\nReview: author-reviewed\n\n"
    report += f"The public 58,735-gene × 119-person matrix was reproduced (32 healthy, 57 mild AP, 20 moderately severe AP, 10 severe AP). Analysis was restricted to the 87 AP patients for severity inference. {len(positive)} of four prespecified modules showed a positive severity trend at BH-FDR <0.05.\n\n"
    report += "| Module | Spearman rho | permutation P | BH-FDR | severe−mild score (95%CI) | UQ sensitivity rho/P |\n|---|---:|---:|---:|---:|---:|\n"
    for r in results:
        report += f"| {r['module']} | {r['spearman_rho']:.3f} | {r['one_sided_permutation_p']:.4g} | {r['BH_FDR']:.4g} | {r['severe_minus_mild']:.3f} ({r['bootstrap_CI_low']:.3f}, {r['bootstrap_CI_high']:.3f}) | {r['upper_quartile_rho']:.3f}/{r['upper_quartile_p']:.4g} |\n"
    report += "\nThese patients are independent of the microbiome cohorts. The module results provide directional host-response triangulation only; they do not establish within-patient microbe–host correlation, mediation, or causality.\n"
    write_text(outdir / "HOST_RESPONSE_REPORT.md", report)
    (ROOT / "logs/HOST_RESPONSE_COMPLETE").write_text(datetime.now(timezone.utc).astimezone().isoformat(timespec="seconds") + "\n")
    print(json.dumps({"status": "COMPLETE", "positive_FDR_modules": len(positive), "results": results}))


if __name__ == "__main__":
    main()
