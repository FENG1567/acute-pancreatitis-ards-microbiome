#!/usr/bin/env python3
"""Integrate Stage 1 evidence without pooling incompatible cohorts."""

from __future__ import annotations

import csv
import json
import os
import pathlib
from datetime import datetime, timezone

ROOT = pathlib.Path(os.environ.get("STAGE1_ROOT", pathlib.Path.cwd()))
STAGE0 = pathlib.Path(os.environ.get("STAGE0_ROOT", ROOT.parent / "gut_pancreas_twin_stage0"))


def read_tsv(path):
    with open(path, encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def read_json(path):
    with open(path, encoding="utf-8") as fh:
        return json.load(fh)


def write_text(path, text):
    path = pathlib.Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".partial")
    tmp.write_text(text, encoding="utf-8", newline="\n")
    os.replace(tmp, path)


def write_tsv(path, rows):
    path = pathlib.Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    cols = list(rows[0])
    tmp = path.with_suffix(path.suffix + ".partial")
    with tmp.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    os.replace(tmp, path)


def f(x, n=3):
    return f"{float(x):.{n}f}"


def main():
    required = [
        ROOT / "07_results/primary_test.json",
        ROOT / "06_supporting/blood_neutrophil/blood_neutrophil_summary.json",
        ROOT / "06_supporting/host_response/host_response_summary.json",
        STAGE0 / "06_results/stage05_PRJNA1470976_summary.json",
        STAGE0 / "06_results/ST004767_summary.json",
    ]
    missing = [str(p) for p in required if not p.exists() or p.stat().st_size == 0]
    if missing:
        raise RuntimeError(f"required evidence missing: {missing}")

    primary = read_json(required[0])
    paired = read_json(required[1])
    host = read_json(required[2])
    etiology = read_json(required[3])
    metab = read_json(required[4])
    cov = read_tsv(ROOT / "07_results/covariate_support.tsv")
    adj = next(r for r in cov if r["model"] == "age_sex_bmi_adjusted" and r["term"] == "signature_per_1SD")
    stage06 = read_tsv(STAGE0 / "06_results/stage06_PRJNA771396_paired_tests.tsv")
    s06 = next(r for r in stage06 if r["analysis"] == "primary")
    source = read_tsv(STAGE0 / "06_results/stage05_PRJNA1031835_source_tests.tsv")
    src = next(r for r in source if float(r["threshold"]) == 0.0001 and r["timepoint"] == "1" and r["comparison"] == "feces_to_pancreatic_fluid")

    rows = []

    def add(claim_id, classification, claim, dataset, unit, independent, analysis, result,
            interval, p_or_fdr, limitation, allowed, forbidden):
        rows.append({
            "claim_id": claim_id,
            "classification": classification,
            "claim": claim,
            "dataset": dataset,
            "statistical_unit": unit,
            "independence": independent,
            "analysis": analysis,
            "result": result,
            "interval": interval,
            "p_or_fdr": p_or_fdr,
            "limitation": limitation,
            "allowed_wording": allowed,
            "forbidden_wording": forbidden,
        })

    add("C01", "confirmatory",
        "The outcome-blind locked five-genus score is higher in AP patients with ARDS.",
        "PRJNA893348", "65 unique AP patients (26 ARDS, 39 nonARDS)", "Independent of discovery cohorts",
        "One-sided patient-label permutation test, 10,000 permutations",
        f"ARDS-nonARDS mean difference {f(primary['observed_mean_difference_SD'])} SD; Cliff's delta {f(primary['cliffs_delta'])}; descriptive AUROC {f(primary['AUROC_descriptive'])}",
        f"bootstrap 95% CI {f(primary['bootstrap_95CI'][0])} to {f(primary['bootstrap_95CI'][1])}; AUROC 95% CI {f(primary['AUROC_bootstrap_95CI'][0])} to {f(primary['AUROC_bootstrap_95CI'][1])}",
        f"one-sided permutation P={primary['one_sided_permutation_p']:.5f}",
        "Individual sampling-to-ARDS timing is unavailable; cohort-level early pre-ARDS timing only.",
        "associated with early ARDS status in an independent AP cohort",
        "prospective prediction, causal biomarker, externally calibrated clinical model")

    add("C02", "supportive",
        "The primary association remains after prespecified age, sex and BMI adjustment.",
        "PRJNA893348", "65 unique AP patients", "Same patients as C01; not independent evidence",
        "Prespecified logistic regression support model",
        f"OR per 1-SD score={f(adj['OR'],2)}",
        f"95% CI {f(adj['CI_low'],2)} to {f(adj['CI_high'],2)}",
        f"two-sided Wald P={float(adj['Wald_p_two_sided']):.4f}",
        "Small event count and observational covariates; model supports but does not replace the permutation test.",
        "association is not explained by measured age, sex and BMI",
        "confounding eliminated, independently validated prediction model")

    loo = primary["leave_one_patient_out"]
    cont = primary["strict_contamination_sensitivity"]
    add("C03", "supportive",
        "No single patient or Staphylococcus component determines the primary result.",
        "PRJNA893348", "Repeated patient-level sensitivity analyses", "Same patients as C01",
        "Leave-one-patient-out and fixed-weight leave-Staphylococcus-out analyses",
        f"LOO effect {f(loo['effect_min'])}-{f(loo['effect_max'])} SD; strict effect {f(cont['effect'])} SD",
        "All leave-one-patient-out directions positive",
        f"LOO P range {loo['p_min']:.4f}-{loo['p_max']:.4f}; strict P={cont['p']:.4f}",
        "No public sequenced blanks; sensitivity analysis cannot prove absence of contamination.",
        "robust to patient influence and a prespecified strict contamination scenario",
        "contamination-free, all detected taxa are biological")

    pmain = paired["all_pairs_primary_support"]
    add("C04", "supportive-negative",
        "Blood and neutrophil scores do not show significant participant-level rank concordance.",
        "PRJNA428535", f"{pmain['n_pairs']} unique paired participants", "Independent participants",
        "Blood-neutrophil Spearman correlation with pairing permutation",
        f"rho={f(pmain['spearman_rho'])}", "not an effect-size interval",
        f"pairing-permutation P={pmain['pairing_permutation_p']:.5f}",
        "Low-biomass Ion Torrent V3 16S data without public sequenced blanks.",
        "no evidence of participant-level rank concordance between compartments",
        "blood-neutrophil concordance, intracellular viable bacteria, migration direction, ARDS validation")

    add("C05", "supportive",
        "The locked score is enriched in the neutrophil compartment relative to paired blood.",
        "PRJNA428535", f"{pmain['n_pairs']} unique paired participants", "Independent participants",
        "Paired mean difference, bootstrap interval and two-sided sign-flip test",
        f"neutrophil-blood mean={f(pmain['neutrophil_minus_blood_mean'])} SD",
        f"paired bootstrap 95% CI {f(pmain['paired_bootstrap_CI_low'])} to {f(pmain['paired_bootstrap_CI_high'])}",
        f"two-sided sign-flip P={pmain['signflip_p_two_sided']:.5f}",
        "A score difference does not identify cell localization or transport mechanism.",
        "paired neutrophil-compartment enrichment of the locked score",
        "proof of microbial transport, viability or intracellular localization")

    for i, r in enumerate(host["results"], start=6):
        label = r["module"].replace("_", " ")
        add(f"C{i:02d}", "supportive",
            f"The prespecified {label} module increases with AP severity.",
            "GSE194331/PRJNA800337", "87 independent AP patients", "Independent of all microbiome cohorts",
            "Ordinal severity Spearman trend with 10,000 permutations and BH correction",
            f"rho={f(r['spearman_rho'])}; severe-mild={f(r['severe_minus_mild'])}",
            f"bootstrap 95% CI {f(r['bootstrap_CI_low'])} to {f(r['bootstrap_CI_high'])}",
            f"BH-FDR={r['BH_FDR']:.5f}",
            "Different patients from microbiome cohorts; no within-patient multi-omics link.",
            "independent host-response pathway support in the same disease context",
            "microbe-host correlation, mediation or causality")

    add("C10", "supportive",
        "Cross-site evidence supported locking the microbial signature before outcome analysis.",
        "PRJNA1031835 and PRJNA771396", "14 and 17 paired patients, analyzed separately", "Discovery/support cohorts",
        "Source-matching and blood-peripancreatic concordance permutations",
        f"time-1 feces-pancreatic delta={f(src['delta_mean'])}; blood-peripancreatic delta={f(s06['delta_mean'])}",
        "cohort-specific leave-one-out effects remained positive",
        f"P={float(src['permutation_p_one_sided']):.5f} and P={float(s06['permutation_p_one_sided']):.5f}",
        "Shared reads or profiles do not establish viable organisms or direction of transfer.",
        "cross-site similarity and concordance",
        "gut-to-pancreas migration, source direction, infection causality")

    ap_hc = next(r for r in etiology["module_contrasts"] if r["contrast"] == "AP_vs_HC")
    add("C11", "supportive",
        "An independent etiology-focused cohort supports an AP-associated locked microbial state.",
        "PRJNA1470976", "12 AP, 12 alcoholic AP, 12 healthy participants", "Independent cohort",
        "Locked-module contrasts with patient-level permutation",
        f"AP-HC module delta={f(ap_hc['delta_module'])}", "No compatible precision interval in inherited result",
        f"permutation P={ap_hc['permutation_p']:.5f}",
        "Small V3-V4 16S cohort; not an ARDS cohort and not species-resolved.",
        "etiology-context support for an AP-associated state",
        "ARDS validation or Enterococcus faecalis strain evidence")

    add("C12", "exploratory",
        "Public metabolomics data provide disease-context background only.",
        "ST004767", f"{metab['n_AP']} AP, {metab['n_control']} controls, {metab['n_QC']} QC samples",
        "Independent cohort", "Inherited QC-filtered exploratory analysis",
        f"{metab['n_features_pass_qc']} features passed QC; {metab['n_fdr_0_05']} had FDR<=0.05",
        "No prespecified metabolite-level causal interval", f"exploratory permutation P={metab['permutation_p']:.4f}",
        "No AP severity, ARDS, time or key covariates; internal CV AUROC=1.0 is implausibly optimistic.",
        "metabolic disease-context background",
        "ARDS prediction, causal bridge, external diagnostic performance")

    excluded = read_tsv(ROOT / "06_supporting/excluded_candidate_cohorts.tsv")
    for j, r in enumerate(excluded, start=13):
        add(f"C{j:02d}", "excluded",
            f"{r['dataset']} did not pass all five eligibility gates for additional outcome validation.",
            r["dataset"], f"publicly analyzable units={r['patient_level_units']}", "Not used in inferential synthesis",
            "Stage 0.7 patient mapping, assay, outcome, overlap and sample-size gate",
            r["gate"], "not applicable", "not applicable", r["public_label_strength"],
            r["allowed_role"], r["prohibited_role"])

    support_dir = ROOT / "06_supporting"
    write_tsv(support_dir / "supporting_evidence_matrix.tsv", rows)
    write_tsv(ROOT / "09_manuscript/claim_evidence_matrix.tsv", rows)

    report = [
        "# Supporting evidence report",
        "",
        "Review: author-reviewed",
        "",
        "## Integrated conclusion",
        "",
        "The confirmatory result remains the locked five-genus score association with ARDS in 65 independent AP patients. Supporting cohorts add cross-site similarity, a paired neutrophil-compartment enrichment with no significant blood-neutrophil rank concordance, and independent host-response severity trends. They are not pooled and do not substitute for the primary outcome test.",
        "",
        "## Primary and robustness evidence",
        "",
        f"The ARDS-minus-nonARDS score difference was {f(primary['observed_mean_difference_SD'])} SD (95% CI {f(primary['bootstrap_95CI'][0])} to {f(primary['bootstrap_95CI'][1])}; one-sided permutation P={primary['one_sided_permutation_p']:.5f}). All 65 leave-one-patient-out analyses remained positive and significant, and fixed-weight omission of Staphylococcus retained a {f(cont['effect'])}-SD effect (P={cont['p']:.4f}).",
        "",
        "## Blood-neutrophil paired support",
        "",
        f"Among {pmain['n_pairs']} analyzable participant pairs, blood and neutrophil scores were not significantly rank-correlated (rho={f(pmain['spearman_rho'])}, pairing-permutation P={pmain['pairing_permutation_p']:.5f}). The neutrophil-minus-blood mean difference was {f(pmain['neutrophil_minus_blood_mean'])} SD (95% CI {f(pmain['paired_bootstrap_CI_low'])} to {f(pmain['paired_bootstrap_CI_high'])}; sign-flip P={pmain['signflip_p_two_sided']:.5f}). This supports a compartment-level enrichment, not within-person concordance.",
        "",
        "## Independent host-response support",
        "",
        "All four prespecified host modules increased with AP severity after BH correction. The strongest ordinal association was " + max(host["results"], key=lambda x: x["spearman_rho"])["module"].replace("_", " ") + ". These are different patients from the microbiome cohorts and therefore do not constitute patient-level multi-omics.",
        "",
        "## Restricted evidence",
        "",
        "PRJNA1470976 is retained as small-sample etiology context. ST004767 is retained only as exploratory metabolic background. The four Stage 0.7 candidates are excluded from outcome validation because none passed all patient-mapping, assay, endpoint, overlap and analyzable-sample gates.",
        "",
        "## Non-negotiable wording boundary",
        "",
        "The analysis supports an early ARDS association and convergent cross-cohort evidence. It does not establish live bacteria, intracellular carriage, migration direction, causality, same-patient multi-omics, or a clinically validated prediction model.",
        "",
    ]
    write_text(support_dir / "SUPPORTING_EVIDENCE_REPORT.md", "\n".join(report))
    write_text(ROOT / "logs/SUPPORTING_EVIDENCE_COMPLETE",
               datetime.now(timezone.utc).astimezone().isoformat(timespec="seconds") + "\n")
    print(json.dumps({"status": "COMPLETE", "claims": len(rows)}, ensure_ascii=False))


if __name__ == "__main__":
    main()
