#!/usr/bin/env python3
"""Run the protocol-locked PRJNA893348 patient-level primary analysis."""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import math
import os
import pathlib

import numpy as np
from scipy import stats


ROOT = pathlib.Path(__file__).resolve().parents[2]
DATA = ROOT / "data"
OUT = ROOT / "results" / "analysis"
SEED = 20260929
B = 10_000
FEATURES = ["Borreliella", "Escherichia", "Staphylococcus", "Enterococcus", "Klebsiella"]
INPUT_MAP = {"Borreliella": None, "Escherichia": "Escherichia/Shigella",
             "Staphylococcus": "Staphylococcus", "Enterococcus": "Enterococcus",
             "Klebsiella": "Klebsiella"}


def sha256(path: pathlib.Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def read_tsv(path: pathlib.Path) -> list[dict]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write_text(path: pathlib.Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".partial")
    tmp.write_text(text, encoding="utf-8", newline="\n")
    os.replace(tmp, path)


def write_json(path: pathlib.Path, obj) -> None:
    write_text(path, json.dumps(obj, ensure_ascii=False, indent=2, allow_nan=False) + "\n")


def clean(v):
    if isinstance(v, np.generic):
        return v.item()
    if isinstance(v, float) and not math.isfinite(v):
        return None
    if isinstance(v, dict):
        return {k: clean(x) for k, x in v.items()}
    if isinstance(v, list):
        return [clean(x) for x in v]
    return v


def write_tsv(path: pathlib.Path, rows: list[dict], columns: list[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if columns is None:
        columns = list(rows[0]) if rows else []
    tmp = path.with_suffix(path.suffix + ".partial")
    with tmp.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=columns, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows([{k: clean(v) for k, v in r.items()} for r in rows])
    os.replace(tmp, path)


def score_from_counts(counts: np.ndarray, genera: list[str], pseudocount: float,
                      omitted: set[str] | None = None) -> tuple[np.ndarray, dict[str, np.ndarray]]:
    omitted = omitted or set()
    g = list(genera)
    x = counts.astype(float)
    if "Borreliella" not in g:
        g.append("Borreliella")
        x = np.column_stack([x, np.zeros(x.shape[0])])
    logs = np.log(x + pseudocount)
    clr = logs - logs.mean(axis=1, keepdims=True)
    components = {}
    total = np.zeros(x.shape[0])
    for f in FEATURES:
        col = INPUT_MAP[f] if INPUT_MAP[f] is not None else "Borreliella"
        arr = clr[:, g.index(col)]
        components[f] = arr
        if f not in omitted:
            total += 0.2 * arr
    sd = total.std(ddof=1)
    if not sd > 0:
        raise RuntimeError("zero variance signature score")
    z = (total - total.mean()) / sd
    return z, components


def auc(y: np.ndarray, s: np.ndarray) -> float:
    pos = s[y == 1]
    neg = s[y == 0]
    return float(np.mean(pos[:, None] > neg[None, :]) + 0.5 * np.mean(pos[:, None] == neg[None, :]))


def effect(y: np.ndarray, s: np.ndarray) -> float:
    return float(s[y == 1].mean() - s[y == 0].mean())


def perm_test(y: np.ndarray, s: np.ndarray, seed: int, b: int = B) -> tuple[float, np.ndarray]:
    rng = np.random.default_rng(seed)
    obs = effect(y, s)
    vals = np.empty(b)
    for i in range(b):
        vals[i] = effect(rng.permutation(y), s)
    p = (1.0 + np.sum(vals >= obs)) / (b + 1.0)
    return float(p), vals


def strat_boot(y: np.ndarray, s: np.ndarray, seed: int, b: int = B) -> tuple[np.ndarray, np.ndarray]:
    rng = np.random.default_rng(seed)
    ip = np.flatnonzero(y == 1)
    ineg = np.flatnonzero(y == 0)
    e = np.empty(b)
    a = np.empty(b)
    for i in range(b):
        bp = rng.choice(ip, len(ip), replace=True)
        bn = rng.choice(ineg, len(ineg), replace=True)
        idx = np.concatenate([bp, bn])
        yy = np.concatenate([np.ones(len(bp), int), np.zeros(len(bn), int)])
        ss = s[idx]
        e[i] = effect(yy, ss)
        a[i] = auc(yy, ss)
    return e, a


def logistic_fit(y: np.ndarray, X: np.ndarray, maxiter: int = 100) -> dict:
    X = np.asarray(X, float)
    y = np.asarray(y, float)
    beta = np.zeros(X.shape[1])
    converged = False
    for it in range(maxiter):
        eta = np.clip(X @ beta, -35, 35)
        p = 1 / (1 + np.exp(-eta))
        w = np.maximum(p * (1 - p), 1e-9)
        info = X.T @ (w[:, None] * X)
        score = X.T @ (y - p)
        ridge = np.eye(X.shape[1]) * 1e-9
        ridge[0, 0] = 0
        try:
            step = np.linalg.solve(info + ridge, score)
        except np.linalg.LinAlgError:
            step = np.linalg.pinv(info + ridge) @ score
        beta += step
        if np.max(np.abs(step)) < 1e-8:
            converged = True
            break
    eta = np.clip(X @ beta, -35, 35)
    p = 1 / (1 + np.exp(-eta))
    w = np.maximum(p * (1 - p), 1e-9)
    cov = np.linalg.pinv(X.T @ (w[:, None] * X))
    se = np.sqrt(np.diag(cov))
    z = beta / se
    pv = 2 * stats.norm.sf(np.abs(z))
    return {"beta": beta, "se": se, "p": pv, "converged": converged,
            "iterations": it + 1, "predicted": p}


def bh(pvals: list[float]) -> list[float]:
    p = np.asarray(pvals, float)
    order = np.argsort(p)
    q = np.empty(len(p))
    running = 1.0
    for rank_index in range(len(p) - 1, -1, -1):
        idx = order[rank_index]
        val = p[idx] * len(p) / (rank_index + 1)
        running = min(running, val)
        q[idx] = min(1.0, running)
    return q.tolist()


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    meta = read_tsv(DATA / "metadata" / "subject_sample_run_map.tsv")
    meta = [r for r in meta if r["include_primary"] == "true"]
    counts_rows = read_tsv(DATA / "processed" / "PRJNA893348_genus_counts.tsv.gz")
    cb = {r["sample_accession"]: r for r in counts_rows}
    genera = [g for g in counts_rows[0] if g != "sample_accession"]
    meta.sort(key=lambda r: r["subject_id"])
    ids = [r["subject_id"] for r in meta]
    samples = [r["sample_id"] for r in meta]
    y = np.asarray([int(r["ARDS_label_original"]) for r in meta], int)
    Xcount = np.asarray([[int(float(cb[s][g] or 0)) for g in genera] for s in samples], dtype=float)
    depth = Xcount.sum(axis=1)
    score, components = score_from_counts(Xcount, genera, 0.5)
    obs = effect(y, score)
    p_main, perms = perm_test(y, score, SEED + 101)
    mcse = math.sqrt(p_main * (1 - p_main) / (B + 1))
    boot_eff, boot_auc = strat_boot(y, score, SEED + 202)
    eff_ci = np.quantile(boot_eff, [0.025, 0.975])
    a = auc(y, score)
    auc_ci = np.quantile(boot_auc, [0.025, 0.975])
    med_diff = float(np.median(score[y == 1]) - np.median(score[y == 0]))
    cliff = float(2 * a - 1)

    component_rows = []
    for i, r in enumerate(meta):
        row = {"patient_id": ids[i], "sample_id": samples[i], "ARDS": int(y[i]),
               "age": float(r["age"]), "sex": r["sex"], "bmi": float(r["bmi"]),
               "genus_assigned_depth": int(depth[i]), "signature_score_sd": float(score[i]),
               "technical_replicate_handling": "none among AP primary patients",
               "missing_flag": "false", "analysis_population": "primary_65_AP"}
        for f in FEATURES:
            row[f"clr_{f}"] = float(components[f][i])
        component_rows.append(row)
    write_tsv(OUT / "patient_signature_components.tsv", component_rows)

    # Full permutation values make the result directly auditable.
    write_tsv(OUT / "permutation_summary.tsv",
              [{"permutation_id": i + 1, "mean_difference": float(v),
                "at_least_observed": int(v >= obs)} for i, v in enumerate(perms)])

    # Leave-one-patient-out with independent deterministic streams.
    loo = []
    for i, pid in enumerate(ids):
        keep = np.arange(len(ids)) != i
        e = effect(y[keep], score[keep])
        pp, _ = perm_test(y[keep], score[keep], SEED + 1000 + i)
        loo.append({"omitted_patient": pid, "omitted_outcome": int(y[i]),
                    "effect_mean_difference": e, "one_sided_permutation_p": pp,
                    "effect_change_from_full": e - obs, "direction_positive": str(e > 0).lower()})
    write_tsv(OUT / "leave_one_patient_out.tsv", loo)

    # Outcome-blind residualisation against log depth, and a blind low-depth exclusion.
    ld = np.log(depth)
    Z = np.column_stack([np.ones(len(ld)), (ld - ld.mean()) / ld.std(ddof=1)])
    residual = score - Z @ np.linalg.lstsq(Z, score, rcond=None)[0]
    residual = (residual - residual.mean()) / residual.std(ddof=1)
    p_depth, _ = perm_test(y, residual, SEED + 303)
    depth_effect = effect(y, residual)
    all_blind = read_tsv(DATA / "metadata" / "blind_sample_qc.tsv")
    q05 = float(np.quantile([float(r["genus_assigned_depth"]) for r in all_blind], 0.05))
    keep_depth = depth >= q05
    p_low, _ = perm_test(y[keep_depth], score[keep_depth], SEED + 304)

    # Predeclared alternate pseudocount.
    score_pc1, _ = score_from_counts(Xcount, genera, 1.0)
    p_pc1, _ = perm_test(y, score_pc1, SEED + 305)

    sensitivity = [
        {"analysis": "primary_pseudocount_0.5", "n": len(y), "effect": obs, "p": p_main, "critical": "true"},
        {"analysis": "log_depth_residualized", "n": len(y), "effect": depth_effect, "p": p_depth, "critical": "true"},
        {"analysis": "exclude_blind_bottom_5pct_depth", "n": int(keep_depth.sum()), "effect": effect(y[keep_depth], score[keep_depth]), "p": p_low, "critical": "true"},
        {"analysis": "pseudocount_1.0", "n": len(y), "effect": effect(y, score_pc1), "p": p_pc1, "critical": "true"},
    ]
    for j, f in enumerate(FEATURES):
        ss, _ = score_from_counts(Xcount, genera, 0.5, omitted={f})
        pp, _ = perm_test(y, ss, SEED + 400 + j)
        sensitivity.append({"analysis": f"leave_out_{f}_fixed_weight", "n": len(y),
                            "effect": effect(y, ss), "p": pp, "critical": "true"})
    write_tsv(OUT / "contamination_sensitivity.tsv", sensitivity)

    # Exploratory single-feature diagnostics; these never replace the total score.
    genus_rows = []
    for j, f in enumerate(FEATURES):
        arr = components[f]
        if arr.std(ddof=1) > 0:
            arr = (arr - arr.mean()) / arr.std(ddof=1)
            pp, _ = perm_test(y, arr, SEED + 500 + j)
            ee = effect(y, arr)
        else:
            pp, ee = 1.0, 0.0
        genus_rows.append({"genus": f, "effect_sd": ee, "one_sided_p": pp})
    qs = bh([r["one_sided_p"] for r in genus_rows])
    for r, q in zip(genus_rows, qs):
        r["BH_FDR"] = q
    write_tsv(OUT / "exploratory_single_genus.tsv", genus_rows)

    # Prespecified supportive logistic models.
    age = np.asarray([float(r["age"]) for r in meta])
    bmi = np.asarray([float(r["bmi"]) for r in meta])
    male = np.asarray([1.0 if r["sex"].lower() == "male" else 0.0 for r in meta])
    agez = (age - age.mean()) / age.std(ddof=1)
    bmiz = (bmi - bmi.mean()) / bmi.std(ddof=1)
    models = []
    for name, XX, terms in [
        ("unadjusted", np.column_stack([np.ones(len(y)), score]), ["intercept", "signature_per_1SD"]),
        ("age_sex_bmi_adjusted", np.column_stack([np.ones(len(y)), score, agez, male, bmiz]),
         ["intercept", "signature_per_1SD", "age_per_1SD", "male_vs_female", "BMI_per_1SD"]),
        ("age_sex_bmi_logdepth_adjusted", np.column_stack([np.ones(len(y)), score, agez, male, bmiz, Z[:, 1]]),
         ["intercept", "signature_per_1SD", "age_per_1SD", "male_vs_female", "BMI_per_1SD", "log_depth_per_1SD"]),
    ]:
        fit = logistic_fit(y, XX)
        for k, term in enumerate(terms):
            models.append({"model": name, "term": term, "beta": float(fit["beta"][k]),
                           "SE": float(fit["se"][k]), "OR": float(np.exp(fit["beta"][k])),
                           "CI_low": float(np.exp(fit["beta"][k] - 1.96 * fit["se"][k])),
                           "CI_high": float(np.exp(fit["beta"][k] + 1.96 * fit["se"][k])),
                           "Wald_p_two_sided": float(fit["p"][k]),
                           "converged": str(fit["converged"]).lower(), "iterations": fit["iterations"]})
    write_tsv(OUT / "covariate_support.tsv", models)

    # Influence diagnostic: leave-one-out standardized change in the primary mean difference.
    se_boot = float(np.std(boot_eff, ddof=1))
    influence = []
    for row in loo:
        influence.append({"patient_id": row["omitted_patient"], "outcome": row["omitted_outcome"],
                          "loo_effect": row["effect_mean_difference"],
                          "delta_effect": row["effect_change_from_full"],
                          "standardized_delta": row["effect_change_from_full"] / se_boot if se_boot else 0.0,
                          "high_influence_abs_gt_1": str(abs(row["effect_change_from_full"] / se_boot) > 1).lower() if se_boot else "false"})
    write_tsv(OUT / "influence_diagnostics.tsv", influence)

    loo_all_sig = all(float(r["one_sided_permutation_p"]) < 0.05 for r in loo)
    loo_all_positive = all(r["direction_positive"] == "true" for r in loo)
    strict = next(r for r in sensitivity if r["analysis"] == "leave_out_Staphylococcus_fixed_weight")

    primary = {
        "n": int(len(y)), "ARDS": int(y.sum()), "nonARDS": int((1-y).sum()),
        "observed_mean_difference_SD": obs, "bootstrap_95CI": eff_ci.tolist(),
        "median_difference_SD": med_diff, "cliffs_delta": cliff,
        "one_sided_permutation_p": p_main, "permutations": B,
        "monte_carlo_SE": mcse, "minimum_attainable_p": 1 / (B + 1),
        "AUROC_descriptive": a, "AUROC_bootstrap_95CI": auc_ci.tolist(),
        "leave_one_patient_out": {
            "effect_min": min(float(r["effect_mean_difference"]) for r in loo),
            "effect_max": max(float(r["effect_mean_difference"]) for r in loo),
            "p_min": min(float(r["one_sided_permutation_p"]) for r in loo),
            "p_max": max(float(r["one_sided_permutation_p"]) for r in loo),
            "all_direction_positive": loo_all_positive, "all_p_lt_0.05": loo_all_sig,
        },
        "strict_contamination_sensitivity": strict,
        "individual_timing": "unavailable; cohort-level early pre-ARDS timing only",
        "source_lock_sha256": json.loads((DATA / "frozen" / "signature_lock.json").read_text(encoding="utf-8"))["source_signature_sha256"],
    }
    write_json(OUT / "primary_test.json", clean(primary))
    write_tsv(OUT / "primary_effects.tsv", [{
        "effect": "ARDS_minus_nonARDS_mean_signature_SD", "estimate": obs,
        "CI_low": eff_ci[0], "CI_high": eff_ci[1], "one_sided_permutation_p": p_main,
        "AUROC_descriptive": a, "AUROC_CI_low": auc_ci[0], "AUROC_CI_high": auc_ci[1]
    }])

    run_inputs = [DATA / "frozen" / "protocol_lock_manifest.tsv",
                  DATA / "processed" / "PRJNA893348_genus_counts.tsv.gz",
                  DATA / "metadata" / "subject_sample_run_map.tsv"]
    run_manifest = {
        "started_after_protocol_lock": True,
        "master_seed": SEED, "n_permutations": B, "patient_ids": ids,
        "outcome_counts": {"ARDS": int(y.sum()), "nonARDS": int((1-y).sum())},
        "inputs": [{"path": p.relative_to(ROOT).as_posix(), "sha256": sha256(p), "bytes": p.stat().st_size} for p in run_inputs],
        "script_sha256": sha256(pathlib.Path(__file__)),
        "git_protocol_commit": "not-applicable-public-reproduction",
    }
    write_json(OUT / "primary_reproduction_manifest.json", run_manifest)

    # Public reproduction ends with the machine-readable statistical outputs
    # required to reproduce the manuscript results.
    print(json.dumps(clean(primary), ensure_ascii=False))


if __name__ == "__main__":
    main()
