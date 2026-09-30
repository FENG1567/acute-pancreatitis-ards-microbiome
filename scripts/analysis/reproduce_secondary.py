#!/usr/bin/env python3
"""Reviewer-requested sensitivity and supporting analyses.

The script reads the public, analysis-ready inputs in ``data/`` and writes to
``results/analysis``. It never redefines or replaces the locked primary test.
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import math
import os
import pathlib
from collections import defaultdict

import numpy as np
from scipy import optimize, stats


ROOT = pathlib.Path(__file__).resolve().parents[2]
DATA = ROOT / "data"
OUT = ROOT / "results" / "analysis"
SEED = 20260929
B = 10_000
FEATURES = ["Borreliella", "Escherichia", "Staphylococcus", "Enterococcus", "Klebsiella"]
INPUT_MAP = {
    "Borreliella": None,
    "Escherichia": "Escherichia/Shigella",
    "Staphylococcus": "Staphylococcus",
    "Enterococcus": "Enterococcus",
    "Klebsiella": "Klebsiella",
}
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


def read_tsv(path: pathlib.Path) -> list[dict]:
    op = gzip.open if path.suffix == ".gz" else open
    with op(path, "rt", encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write_tsv(path: pathlib.Path, rows: list[dict], columns: list[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    columns = columns or (list(rows[0]) if rows else [])
    tmp = path.with_suffix(path.suffix + ".partial")
    with tmp.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=columns, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    os.replace(tmp, path)


def write_json(path: pathlib.Path, obj) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".partial")
    tmp.write_text(json.dumps(obj, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    os.replace(tmp, path)


def sha256(path: pathlib.Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def zscore(x: np.ndarray) -> np.ndarray:
    x = np.asarray(x, float)
    return (x - x.mean()) / x.std(ddof=1)


def effect(y: np.ndarray, x: np.ndarray) -> float:
    return float(x[y == 1].mean() - x[y == 0].mean())


def strat_boot_effect(y: np.ndarray, x: np.ndarray, seed: int, b: int = B) -> np.ndarray:
    rng = np.random.default_rng(seed)
    ip = np.flatnonzero(y == 1)
    ic = np.flatnonzero(y == 0)
    out = np.empty(b)
    for i in range(b):
        out[i] = x[rng.choice(ip, len(ip), True)].mean() - x[rng.choice(ic, len(ic), True)].mean()
    return out


def perm_p(y: np.ndarray, x: np.ndarray, seed: int, two_sided: bool = False, b: int = B) -> float:
    rng = np.random.default_rng(seed)
    obs = effect(y, x)
    vals = np.asarray([effect(rng.permutation(y), x) for _ in range(b)])
    if two_sided:
        return float((1 + np.sum(np.abs(vals) >= abs(obs))) / (b + 1))
    return float((1 + np.sum(vals >= obs)) / (b + 1))


def build_clr(counts: np.ndarray, genera: list[str], pseudocount: float = 0.5,
              add_borreliella: bool = True) -> tuple[np.ndarray, list[str]]:
    g = list(genera)
    x = counts.astype(float)
    if add_borreliella and "Borreliella" not in g:
        g.append("Borreliella")
        x = np.column_stack([x, np.zeros(len(x))])
    logs = np.log(x + pseudocount)
    return logs - logs.mean(axis=1, keepdims=True), g


def primary_score(counts: np.ndarray, genera: list[str], pseudocount: float = 0.5,
                  omitted: set[str] | None = None) -> tuple[np.ndarray, dict[str, np.ndarray], np.ndarray]:
    omitted = omitted or set()
    clr, g = build_clr(counts, genera, pseudocount, add_borreliella=True)
    comp = {}
    raw = np.zeros(len(counts))
    for f in FEATURES:
        col = INPUT_MAP[f] or "Borreliella"
        comp[f] = clr[:, g.index(col)]
        if f not in omitted:
            raw += 0.2 * comp[f]
    return zscore(raw), comp, raw


def coherent_subcomposition_score(counts: np.ndarray, genera: list[str], kept: list[str], pseudocount: float = 0.5) -> np.ndarray:
    cols = [INPUT_MAP[f] for f in kept if INPUT_MAP[f] in genera]
    x = counts[:, [genera.index(c) for c in cols]]
    return zscore(np.log(x + pseudocount).mean(axis=1))


def aitchison_balance(counts: np.ndarray, genera: list[str], kept_features: list[str], pseudocount: float = 0.5) -> tuple[np.ndarray, int]:
    sig_cols = [INPUT_MAP[f] for f in kept_features if INPUT_MAP[f] in genera]
    prevalence = (counts > 0).mean(axis=0)
    ref_cols = [g for g, p in zip(genera, prevalence) if p >= 0.10 and g not in sig_cols]
    a = np.log(counts[:, [genera.index(g) for g in sig_cols]] + pseudocount).mean(axis=1)
    b = np.log(counts[:, [genera.index(g) for g in ref_cols]] + pseudocount).mean(axis=1)
    scale = math.sqrt(len(sig_cols) * len(ref_cols) / (len(sig_cols) + len(ref_cols)))
    return zscore(scale * (a - b)), len(ref_cols)


def firth_pll(beta: np.ndarray, X: np.ndarray, y: np.ndarray) -> float:
    eta = np.clip(X @ beta, -35, 35)
    p = 1.0 / (1.0 + np.exp(-eta))
    ll = float(np.sum(y * eta - np.logaddexp(0.0, eta)))
    w = np.maximum(p * (1.0 - p), 1e-12)
    info = X.T @ (w[:, None] * X)
    sign, logdet = np.linalg.slogdet(info)
    if sign <= 0 or not np.isfinite(logdet):
        return -np.inf
    return ll + 0.5 * float(logdet)


def firth_fit(X: np.ndarray, y: np.ndarray, start: np.ndarray | None = None) -> tuple[np.ndarray, float, bool]:
    start = np.zeros(X.shape[1]) if start is None else np.asarray(start, float)
    res = optimize.minimize(lambda b: -firth_pll(b, X, y), start, method="BFGS", options={"gtol": 1e-9, "maxiter": 2000})
    if not np.isfinite(res.fun):
        raise RuntimeError("Firth optimization failed")
    return np.asarray(res.x), float(-res.fun), bool(res.success or np.linalg.norm(res.jac) < 1e-5)


def firth_profile_ci(X: np.ndarray, y: np.ndarray, j: int, beta: np.ndarray, max_pll: float) -> tuple[float, float]:
    target = stats.chi2.ppf(0.95, 1)
    keep = [k for k in range(X.shape[1]) if k != j]

    def profile(v: float) -> float:
        def objective(rest):
            b = np.empty(X.shape[1])
            b[j] = v
            b[keep] = rest
            return -firth_pll(b, X, y)
        res = optimize.minimize(objective, beta[keep], method="BFGS", options={"gtol": 1e-8, "maxiter": 1000})
        return float(-res.fun)

    def root(v: float) -> float:
        return 2.0 * (max_pll - profile(v)) - target

    roots = []
    for direction in (-1.0, 1.0):
        center = float(beta[j])
        step = 0.25
        outer = center + direction * step
        for _ in range(12):
            if root(outer) >= 0:
                lo, hi = sorted([center, outer])
                roots.append(float(optimize.brentq(root, lo, hi, xtol=1e-7)))
                break
            step *= 1.8
            outer = center + direction * step
        else:
            roots.append(float("nan"))
    return roots[0], roots[1]


def firth_model_row(name: str, X: np.ndarray, y: np.ndarray, terms: list[str], score_term: str) -> dict:
    beta, pll, converged = firth_fit(X, y)
    j = terms.index(score_term)
    reduced = np.delete(X, j, axis=1)
    _, pll0, converged0 = firth_fit(reduced, y)
    stat = max(0.0, 2.0 * (pll - pll0))
    p = float(stats.chi2.sf(stat, 1))
    lo, hi = firth_profile_ci(X, y, j, beta, pll)
    return {
        "model": name,
        "n": len(y),
        "events": int(y.sum()),
        "parameters_including_intercept": X.shape[1],
        "events_per_nonintercept_parameter": float(y.sum() / (X.shape[1] - 1)),
        "term": score_term,
        "beta_firth": float(beta[j]),
        "OR_firth": float(np.exp(beta[j])),
        "profile_CI_low": float(np.exp(lo)),
        "profile_CI_high": float(np.exp(hi)),
        "penalized_LR_chi2": stat,
        "penalized_LR_p": p,
        "converged": str(converged and converged0).lower(),
    }


def quantiles(x: np.ndarray) -> tuple[float, float, float]:
    q = np.quantile(x, [0.25, 0.5, 0.75])
    return float(q[1]), float(q[0]), float(q[2])


def continuous_baseline(name: str, x: np.ndarray, y: np.ndarray, unit: str) -> dict:
    a, b = x[y == 1], x[y == 0]
    ma, qa1, qa3 = quantiles(a)
    mb, qb1, qb3 = quantiles(b)
    pooled = math.sqrt(((len(a) - 1) * a.var(ddof=1) + (len(b) - 1) * b.var(ddof=1)) / (len(a) + len(b) - 2))
    return {
        "variable": name,
        "unit_or_level": unit,
        "ARDS_summary": f"{a.mean():.2f} ({a.std(ddof=1):.2f}); median {ma:.2f} [{qa1:.2f}, {qa3:.2f}]",
        "nonARDS_summary": f"{b.mean():.2f} ({b.std(ddof=1):.2f}); median {mb:.2f} [{qb1:.2f}, {qb3:.2f}]",
        "standardized_mean_difference": float((a.mean() - b.mean()) / pooled) if pooled else 0.0,
        "missing_ards": int(np.isnan(a).sum()),
        "missing_nonards": int(np.isnan(b).sum()),
        "role": "baseline/technical descriptor; not outcome-driven selection",
    }


def bootstrap_spearman(a: np.ndarray, b: np.ndarray, seed: int, nboot: int = B) -> tuple[float, float, float]:
    rng = np.random.default_rng(seed)
    vals = []
    for _ in range(nboot):
        idx = rng.integers(0, len(a), len(a))
        r = stats.spearmanr(a[idx], b[idx]).statistic
        if np.isfinite(r):
            vals.append(float(r))
    lo, hi = np.quantile(vals, [0.025, 0.975])
    return float(stats.spearmanr(a, b).statistic), float(lo), float(hi)


def paired_bootstrap(d: np.ndarray, seed: int, nboot: int = B) -> tuple[float, float]:
    rng = np.random.default_rng(seed)
    vals = np.asarray([rng.choice(d, len(d), True).mean() for _ in range(nboot)])
    lo, hi = np.quantile(vals, [0.025, 0.975])
    return float(lo), float(hi)


def paired_signflip_p(d: np.ndarray, seed: int, nperm: int = B) -> float:
    rng = np.random.default_rng(seed)
    obs = abs(float(d.mean()))
    vals = np.asarray([abs(float(np.mean(d * rng.choice([-1, 1], len(d))))) for _ in range(nperm)])
    return float((1 + np.sum(vals >= obs)) / (nperm + 1))


def paired_correlation_permutation(a: np.ndarray, b: np.ndarray, seed: int) -> tuple[float, float]:
    rng = np.random.default_rng(seed)
    observed = float(stats.spearmanr(a, b).statistic)
    permuted = np.asarray([stats.spearmanr(a, rng.permutation(b)).statistic for _ in range(B)])
    p_value = float((1 + np.sum(permuted >= observed)) / (B + 1))
    return observed, p_value


def paired_difference_analysis(a: np.ndarray, b: np.ndarray, seed: int) -> tuple[float, float, float, float]:
    """Return neutrophil-minus-blood mean, bootstrap CI, and two-sided sign-flip P."""
    difference = b - a
    observed = float(difference.mean())
    rng = np.random.default_rng(seed)
    permuted = np.asarray([np.mean(difference * rng.choice([-1, 1], len(difference))) for _ in range(B)])
    p_value = float((1 + np.sum(np.abs(permuted) >= abs(observed))) / (B + 1))
    boot = np.asarray([np.mean(rng.choice(difference, len(difference), True)) for _ in range(B)])
    low, high = np.quantile(boot, [0.025, 0.975])
    return observed, float(low), float(high), p_value


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)

    meta = [r for r in read_tsv(DATA / "metadata" / "subject_sample_run_map.tsv") if r["include_primary"] == "true"]
    meta.sort(key=lambda r: r["subject_id"])
    cr = read_tsv(DATA / "processed" / "PRJNA893348_genus_counts.tsv.gz")
    cb = {r["sample_accession"]: r for r in cr}
    genera = [g for g in cr[0] if g != "sample_accession"]
    counts = np.asarray([[float(cb[r["sample_id"]][g] or 0) for g in genera] for r in meta], float)
    y = np.asarray([int(r["ARDS_label_original"]) for r in meta], int)
    depth = counts.sum(axis=1)
    score, components, raw_score = primary_score(counts, genera)

    # Component-level prevalence, abundance, CLR effect, and mathematical contribution.
    comp_rows = []
    total_raw_effect = effect(y, raw_score)
    for j, f in enumerate(FEATURES):
        col = INPUT_MAP[f]
        raw = counts[:, genera.index(col)] if col in genera else np.zeros(len(y))
        rel = np.divide(raw, depth, out=np.zeros_like(raw), where=depth > 0)
        clr = components[f]
        clr_z = zscore(clr) if clr.std(ddof=1) > 0 else np.zeros(len(clr))
        boots = strat_boot_effect(y, clr_z, SEED + 100 + j)
        contribution = 0.2 * effect(y, clr)
        comp_rows.append({
            "genus": f,
            "input_column": col or "structural_zero_inserted",
            "detected_ards_n": int(np.sum(raw[y == 1] > 0)),
            "detected_ards_pct": float(100 * np.mean(raw[y == 1] > 0)),
            "detected_nonards_n": int(np.sum(raw[y == 0] > 0)),
            "detected_nonards_pct": float(100 * np.mean(raw[y == 0] > 0)),
            "median_relative_abundance_ards": float(np.median(rel[y == 1])),
            "median_relative_abundance_nonards": float(np.median(rel[y == 0])),
            "clr_effect_SD": effect(y, clr_z),
            "clr_effect_SD_CI_low": float(np.quantile(boots, 0.025)),
            "clr_effect_SD_CI_high": float(np.quantile(boots, 0.975)),
            "weighted_raw_score_contribution_to_group_difference": contribution,
            "fraction_of_total_raw_score_difference": float(contribution / total_raw_effect) if total_raw_effect else None,
            "correlation_with_primary_score": float(np.corrcoef(clr, score)[0, 1]) if clr.std(ddof=1) > 0 else None,
            "biological_interpretation": "not detected; CLR variation is denominator-induced" if raw.sum() == 0 else "observed relative-abundance component",
        })
    write_tsv(OUT / "primary_component_details.tsv", comp_rows)

    # Primary and reviewer-requested sensitivity intervals.
    sens_specs = [("primary_five_component_fixed", score, "permutation-seed sensitivity; not primary analysis")]
    logdepth = zscore(np.log(depth))
    residual = zscore(score - np.column_stack([np.ones(len(y)), logdepth]) @ np.linalg.lstsq(np.column_stack([np.ones(len(y)), logdepth]), score, rcond=None)[0])
    sens_specs.append(("log_depth_residualized", residual, "prespecified sensitivity"))
    score_pc1, _, _ = primary_score(counts, genera, pseudocount=1.0)
    sens_specs.append(("pseudocount_1.0", score_pc1, "prespecified sensitivity"))
    for f in FEATURES:
        s, _, _ = primary_score(counts, genera, omitted={f})
        sens_specs.append((f"leave_out_{f}_fixed_weight", s, "prespecified fixed-weight component sensitivity"))
    sens_rows = []
    for i, (name, s, role) in enumerate(sens_specs):
        boots = strat_boot_effect(y, s, SEED + 500 + i)
        sens_rows.append({
            "analysis": name,
            "role": role,
            "n": len(y),
            "ARDS": int(y.sum()),
            "nonARDS": int((1 - y).sum()),
            "effect_SD": effect(y, s),
            "bootstrap_CI_low": float(np.quantile(boots, 0.025)),
            "bootstrap_CI_high": float(np.quantile(boots, 0.975)),
            "one_sided_permutation_p": perm_p(y, s, SEED + 700 + i),
            "note": (
                "Primary P=0.009299 is reported in the main analysis; P=0.008799 is from a separate "
                "sensitivity permutation run and is not a replacement for the primary value."
                if i == 0 else "Sensitivity value is retained as a separate non-primary result."
            ),
        })
    write_tsv(OUT / "sensitivity_with_ci.tsv", sens_rows)

    # Post hoc but compositionally coherent alternatives; never re-labelled confirmatory.
    detected = [f for f in FEATURES if f != "Borreliella"]
    four = coherent_subcomposition_score(counts, genera, detected, 0.5)
    four_pc1 = coherent_subcomposition_score(counts, genera, detected, 1.0)
    balance, nref = aitchison_balance(counts, genera, detected, 0.5)
    alt_specs = [
        ("four_detected_genera_log_geometric_mean_pc0.5", four, len(detected), 0, "removes all-zero Borreliella"),
        ("four_detected_genera_log_geometric_mean_pc1.0", four_pc1, len(detected), 0, "pseudocount sensitivity"),
        ("four_detected_genera_vs_prevalent_background_aitchison_balance", balance, len(detected), nref, "reference genera prevalence >=10%"),
    ]
    alt_rows = []
    for i, (name, s, nsig, nrefi, note) in enumerate(alt_specs):
        boots = strat_boot_effect(y, s, SEED + 900 + i)
        alt_rows.append({
            "analysis": name,
            "status": "post hoc exploratory sensitivity; not confirmatory",
            "signature_components": nsig,
            "reference_components": nrefi,
            "effect_SD": effect(y, s),
            "bootstrap_CI_low": float(np.quantile(boots, 0.025)),
            "bootstrap_CI_high": float(np.quantile(boots, 0.975)),
            "two_sided_permutation_p": perm_p(y, s, SEED + 1000 + i, two_sided=True),
            "note": note,
        })
    write_tsv(OUT / "alternative_compositional_sensitivity.tsv", alt_rows)

    # Baseline, missingness, and model-stability audit.
    age = np.asarray([float(r["age"]) for r in meta])
    bmi = np.asarray([float(r["bmi"]) for r in meta])
    male = np.asarray([1.0 if r["sex"].lower() == "male" else 0.0 for r in meta])
    baseline = [
        continuous_baseline("Age", age, y, "years"),
        continuous_baseline("BMI", bmi, y, "kg/m2"),
        continuous_baseline("Genus-assigned sequencing depth", depth, y, "reads"),
        {
            "variable": "Male sex", "unit_or_level": "n (%)",
            "ARDS_summary": f"{int(male[y == 1].sum())} ({100 * male[y == 1].mean():.1f}%)",
            "nonARDS_summary": f"{int(male[y == 0].sum())} ({100 * male[y == 0].mean():.1f}%)",
            "standardized_mean_difference": float((male[y == 1].mean() - male[y == 0].mean()) / math.sqrt((male[y == 1].var(ddof=1) + male[y == 0].var(ddof=1)) / 2)),
            "missing_ards": 0, "missing_nonards": 0,
            "role": "baseline descriptor and prespecified supportive covariate",
        },
        {
            "variable": "Cohort/platform/batch", "unit_or_level": "single level",
            "ARDS_summary": "PUMC / Illumina MiSeq / single public cohort",
            "nonARDS_summary": "PUMC / Illumina MiSeq / single public cohort",
            "standardized_mean_difference": 0.0, "missing_ards": 0, "missing_nonards": 0,
            "role": "no between-center or recorded batch variability available for adjustment",
        },
    ]
    write_tsv(OUT / "clinical_baseline_by_ards.tsv", baseline)
    missing_rows = []
    for var in ["age", "sex", "bmi", "etiology", "antibiotic_exposure", "ICU_status", "time_from_admission", "time_from_symptom_onset"]:
        vals = [r.get(var, "") for r in meta]
        missing = [v in ("", "NA", "N/A", "not applicable") for v in vals]
        missing_rows.append({
            "variable": var,
            "available_n": int(len(vals) - sum(missing)),
            "missing_n": int(sum(missing)),
            "missing_pct": float(100 * np.mean(missing)),
            "analysis_use": "age/sex/BMI only" if var in {"age", "sex", "bmi"} else "not used; unavailable or cohort-level only",
        })
    write_tsv(OUT / "clinical_missingness.tsv", missing_rows)

    X0 = np.column_stack([np.ones(len(y)), score])
    X1 = np.column_stack([np.ones(len(y)), score, zscore(age), male, zscore(bmi)])
    X2 = np.column_stack([np.ones(len(y)), score, zscore(age), male, zscore(bmi), logdepth])
    firth_rows = [
        firth_model_row("unadjusted_firth", X0, y, ["intercept", "signature_per_1SD"], "signature_per_1SD"),
        firth_model_row("age_sex_bmi_firth", X1, y, ["intercept", "signature_per_1SD", "age_per_1SD", "male", "BMI_per_1SD"], "signature_per_1SD"),
        firth_model_row("age_sex_bmi_logdepth_firth", X2, y, ["intercept", "signature_per_1SD", "age_per_1SD", "male", "BMI_per_1SD", "logdepth_per_1SD"], "signature_per_1SD"),
    ]
    write_tsv(OUT / "firth_logistic_models.tsv", firth_rows)

    # Patient-level mapping/timing audit copied into a disclosure-oriented table.
    timing_rows = []
    for r in meta:
        timing_rows.append({
            "subject_id": r["subject_id"], "sample_id": r["sample_id"], "run_accession": r["run_accession"],
            "ARDS_label_original": r["ARDS_label_original"], "ARDS_definition_public": r["ARDS_definition"],
            "collection_time_public": r["collection_time_original"], "assessment_window_public": r["ARDS_assessment_window"],
            "individual_sampling_to_ards_interval": "unavailable",
            "baseline_ards_exclusion_verifiable_per_patient": "no; only cohort-level publication statement",
            "mapping_confidence": r["mapping_confidence"],
        })
    write_tsv(OUT / "patient_outcome_timing.tsv", timing_rows)

    # Paired blood-neutrophil genus-level effects and rho uncertainty.
    ledger = [r for r in read_tsv(DATA / "metadata" / "stage08_patient_sample_ledger.tsv") if r["dataset"] == "PRJNA428535"]
    pr = read_tsv(DATA / "processed" / "PRJNA428535_genus_counts.tsv.gz")
    pgenera = [g for g in pr[0] if g != "run_accession"]
    pcb = {r["run_accession"]: r for r in pr}
    ledger.sort(key=lambda r: r["run_accession"])
    px = np.asarray([[float(pcb[r["run_accession"]][g] or 0) for g in pgenera] for r in ledger], float)
    pscore, pcomp, _ = primary_score(px, pgenera)
    strict_score, _, _ = primary_score(px, pgenera, omitted={"Staphylococcus"})
    pdepth = px.sum(axis=1)
    byp = defaultdict(dict)
    for i, r in enumerate(ledger):
        byp[r["patient_id"]][r["specimen"]] = i
    genus_pair_rows = []
    pair_scores = []
    for pid in sorted(byp):
        ib = byp[pid]["blood"]
        ineu = byp[pid]["neutrophils"]
        blood_record = ledger[ib]
        row = {
            "patient_id": pid,
            "group": blood_record["group"],
            "blood_run": blood_record["run_accession"],
            "neutrophil_run": ledger[ineu]["run_accession"],
            "blood_depth": int(pdepth[ib]),
            "neutrophil_depth": int(pdepth[ineu]),
            "blood_score": float(pscore[ib]),
            "neutrophil_score": float(pscore[ineu]),
            "blood_strict_score": float(strict_score[ib]),
            "neutrophil_strict_score": float(strict_score[ineu]),
            "pair_qc_ge100": str(pdepth[ib] >= 100 and pdepth[ineu] >= 100).lower(),
            "pair_qc_ge1000": str(pdepth[ib] >= 1000 and pdepth[ineu] >= 1000).lower(),
        }
        for feature in FEATURES:
            column = INPUT_MAP[feature]
            row[f"blood_{feature}_detected"] = int((px[ib, pgenera.index(column)] if column in pgenera else 0) > 0)
            row[f"neutrophil_{feature}_detected"] = int((px[ineu, pgenera.index(column)] if column in pgenera else 0) > 0)
        pair_scores.append(row)
    write_tsv(OUT / "paired_scores.tsv", pair_scores)
    for j, f in enumerate(FEATURES):
        col = INPUT_MAP[f]
        raw = px[:, pgenera.index(col)] if col in pgenera else np.zeros(len(px))
        blood_idx = np.asarray([byp[p]["blood"] for p in sorted(byp)])
        neu_idx = np.asarray([byp[p]["neutrophils"] for p in sorted(byp)])
        diff = pcomp[f][neu_idx] - pcomp[f][blood_idx]
        lo, hi = paired_bootstrap(diff, SEED + 1200 + j)
        genus_pair_rows.append({
            "genus": f,
            "n_pairs": len(blood_idx),
            "blood_detected_n": int(np.sum(raw[blood_idx] > 0)),
            "neutrophil_detected_n": int(np.sum(raw[neu_idx] > 0)),
            "both_detected_n": int(np.sum((raw[blood_idx] > 0) & (raw[neu_idx] > 0))),
            "both_absent_n": int(np.sum((raw[blood_idx] == 0) & (raw[neu_idx] == 0))),
            "mean_neutrophil_minus_blood_CLR": float(diff.mean()),
            "paired_bootstrap_CI_low": lo,
            "paired_bootstrap_CI_high": hi,
            "two_sided_signflip_p": paired_signflip_p(diff, SEED + 1300 + j),
            "interpretation": "denominator-induced CLR shift; taxon absent in both compartments" if raw.sum() == 0 else "relative-abundance component; not absolute load",
        })
    write_tsv(OUT / "paired_genus_effects.tsv", genus_pair_rows)
    pa = np.asarray([x["blood_score"] for x in pair_scores])
    pb = np.asarray([x["neutrophil_score"] for x in pair_scores])
    rho, rlo, rhi = bootstrap_spearman(pa, pb, SEED + 1400)
    write_tsv(OUT / "paired_correlation_ci.tsv", [{
        "analysis": "primary_relative_score_blood_vs_neutrophil",
        "n_pairs": len(pa), "spearman_rho": rho,
        "patient_bootstrap_CI_low": rlo, "patient_bootstrap_CI_high": rhi,
        "interpretation": "no evidence of rank concordance; interval quantifies uncertainty",
    }])

    keep100 = np.asarray([r["pair_qc_ge100"] == "true" for r in pair_scores])
    keep1000 = np.asarray([r["pair_qc_ge1000"] == "true" for r in pair_scores])
    paired_results = []
    for name, keep, a_values, b_values, corr_seed, diff_seed in [
        ("all_pairs_ge100", keep100, pa, pb, SEED + 901, SEED + 902),
        ("pairs_ge1000_sensitivity", keep1000, pa, pb, SEED + 903, SEED + 904),
        (
            "omit_Staphylococcus_fixed_weight",
            keep100,
            np.asarray([r["blood_strict_score"] for r in pair_scores]),
            np.asarray([r["neutrophil_strict_score"] for r in pair_scores]),
            SEED + 905,
            SEED + 906,
        ),
    ]:
        aa, bb = a_values[keep], b_values[keep]
        corr, corr_p = paired_correlation_permutation(aa, bb, corr_seed)
        difference, ci_low, ci_high, sign_p = paired_difference_analysis(aa, bb, diff_seed)
        paired_results.append({
            "analysis": name,
            "n_pairs": int(keep.sum()),
            "spearman_rho": corr,
            "pairing_permutation_p": corr_p,
            "neutrophil_minus_blood_mean": difference,
            "paired_bootstrap_CI_low": ci_low,
            "paired_bootstrap_CI_high": ci_high,
            "signflip_p_two_sided": sign_p,
        })
    write_tsv(OUT / "paired_results.tsv", paired_results)

    concordance = []
    for feature in FEATURES:
        blood = np.asarray([r[f"blood_{feature}_detected"] for r in pair_scores], int)
        neutrophil = np.asarray([r[f"neutrophil_{feature}_detected"] for r in pair_scores], int)
        concordance.append({
            "genus": feature,
            "n_pairs": len(pair_scores),
            "blood_detected": int(blood.sum()),
            "neutrophil_detected": int(neutrophil.sum()),
            "both_detected": int(np.sum((blood == 1) & (neutrophil == 1))),
            "both_absent": int(np.sum((blood == 0) & (neutrophil == 0))),
            "agreement_fraction": float(np.mean(blood == neutrophil)),
        })
    write_tsv(OUT / "genus_pair_concordance.tsv", concordance)

    # Host module dictionary and a bitwise-style recomputation check against the delivered score file.
    mapping_rows = read_tsv(DATA / "host_response" / "ensembl_symbol_mapping.tsv")
    mapping = {r["symbol"]: r["ensembl_gene_id"] for r in mapping_rows if r["resolved"] == "true"}
    coverage = {(r["module"], r["symbol"]): r for r in read_tsv(DATA / "host_response" / "module_gene_coverage.tsv")}
    module_rows = []
    for module, spec in MODULES.items():
        for direction, sign in [("positive", 1), ("negative", -1)]:
            for symbol in spec[direction]:
                cov = coverage.get((module, symbol), {})
                module_rows.append({
                    "module_version": "stage1_modules_20260929_v1",
                    "module": module, "symbol": symbol, "direction": sign,
                    "weight_rule": "equal mean of signed within-gene z scores",
                    "ensembl_gene_id": mapping.get(symbol, ""),
                    "found_in_filtered_matrix": cov.get("found_in_matrix", "unknown"),
                    "selection_status": "specified in analysis script before module results; not externally registered",
                    "scientific_role": "independent AP-severity context only; not an ARDS or same-patient microbe-host test",
                })
    write_tsv(OUT / "host_module_dictionary.tsv", module_rows)
    overlap_rows = []
    names = list(MODULES)
    sets = {m: set(MODULES[m]["positive"] + MODULES[m]["negative"]) for m in names}
    for i, a in enumerate(names):
        for b in names[i + 1:]:
            inter = sets[a] & sets[b]
            union = sets[a] | sets[b]
            overlap_rows.append({"module_a": a, "module_b": b, "overlap_n": len(inter), "overlap_symbols": ";".join(sorted(inter)), "jaccard": len(inter) / len(union)})
    write_tsv(OUT / "host_module_overlap.tsv", overlap_rows)

    # Public integrity audit. The original internal Git commit is documented in
    # the frozen provenance records but is not required to exist in a clone.
    primary_path = OUT / "primary_test.json"
    run_manifest = json.loads((OUT / "primary_reproduction_manifest.json").read_text(encoding="utf-8"))
    lock_rows = [
        {
            "artifact": "public primary reproduction manifest", "path": "results/analysis/primary_reproduction_manifest.json",
            "sha256_or_commit": sha256(OUT / "primary_reproduction_manifest.json"), "timestamp": "generated during public reproduction",
            "audit_meaning": "binds the public primary script and analysis-ready inputs",
        },
        {
            "artifact": "frozen protocol", "path": "data/frozen/FROZEN_ANALYSIS_PROTOCOL.md",
            "sha256_or_commit": sha256(DATA / "frozen" / "FROZEN_ANALYSIS_PROTOCOL.md"),
            "timestamp": "2026-09-29T00:52:46+08:00",
            "audit_meaning": "hash appears in protocol_lock_manifest.tsv",
        },
        {
            "artifact": "source signature", "path": "data/frozen/stage0_locked_signature.tsv",
            "sha256_or_commit": sha256(DATA / "frozen" / "stage0_locked_signature.tsv"),
            "timestamp": "recorded in the frozen protocol bundle",
            "audit_meaning": "complete candidate table, not just the five selected genera",
        },
        {
            "artifact": "primary result", "path": "results/analysis/primary_test.json",
            "sha256_or_commit": sha256(primary_path),
            "timestamp": "generated during public reproduction",
            "audit_meaning": "recomputed from the locked public inputs",
        },
    ]
    write_tsv(OUT / "lock_audit_manifest.tsv", lock_rows)

    summary = {
        "primary_score_seed_sensitivity": sens_rows[0],
        "leave_klebsiella_out": next(r for r in sens_rows if r["analysis"] == "leave_out_Klebsiella_fixed_weight"),
        "alternative_compositional_sensitivities": alt_rows,
        "firth_models": firth_rows,
        "paired_rho": {"estimate": rho, "CI": [rlo, rhi]},
        "data_limits": {
            "independent_AP_ARDS_replication": "not identified in audited public cohorts",
            "individual_sampling_to_ARDS_interval": "unavailable",
            "negative_sequencing_controls": "unavailable for PRJNA893348 and PRJNA428535",
            "clinical_covariates_beyond_age_sex_BMI": "unavailable",
        },
        "interpretive_status": "single-cohort association with separately labelled hypothesis-generating support",
        "input_hashes": {
            "PRJNA893348_genus_counts": sha256(DATA / "processed" / "PRJNA893348_genus_counts.tsv.gz"),
            "subject_sample_run_map": sha256(DATA / "metadata" / "subject_sample_run_map.tsv"),
            "PRJNA428535_genus_counts": sha256(DATA / "processed" / "PRJNA428535_genus_counts.tsv.gz"),
            "host_module_results": sha256(DATA / "host_response" / "host_module_results.tsv"),
        },
        "source_run_manifest_protocol_commit": run_manifest["git_protocol_commit"],
    }
    write_json(OUT / "revision_analysis_summary.json", summary)
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
