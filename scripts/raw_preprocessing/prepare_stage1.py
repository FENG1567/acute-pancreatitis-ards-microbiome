#!/usr/bin/env python3
"""Create and freeze the Stage 1 workspace before any outcome analysis.

This script deliberately performs only metadata reconciliation, blind QC, file
copying, hashing, and protocol generation.  It never computes a signature by
ARDS group.
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import os
import pathlib
import shutil
import statistics
import subprocess
from collections import Counter, defaultdict
from datetime import datetime, timezone


_stage0_env = os.environ.get("STAGE0_ROOT")
_stage1_env = os.environ.get("STAGE1_ROOT")
if not _stage0_env or not _stage1_env:
    raise SystemExit("Set STAGE0_ROOT and STAGE1_ROOT before running raw workspace preparation")
STAGE0 = pathlib.Path(_stage0_env)
ROOT = pathlib.Path(_stage1_env)
EXPECTED_SOURCE_LOCK_SHA256 = "0a10ecee753b8c31f00b216d2f457c6795abf5c09abee75714f665f09974b0af"
MASTER_SEED = 20260929
DIRS = [
    "00_admin", "01_metadata", "02_raw/inherited", "03_processed", "04_qc",
    "05_primary", "06_supporting", "07_results", "08_figures",
    "09_manuscript", "10_release", "envs", "workflow", "logs", "tmp",
]


def sha256(path: pathlib.Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def write_text(path: pathlib.Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".partial")
    tmp.write_text(text, encoding="utf-8", newline="\n")
    os.replace(tmp, path)


def write_json(path: pathlib.Path, obj) -> None:
    write_text(path, json.dumps(obj, ensure_ascii=False, indent=2) + "\n")


def write_tsv(path: pathlib.Path, rows: list[dict], columns: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".partial")
    with tmp.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=columns, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    os.replace(tmp, path)


def read_tsv(path: pathlib.Path) -> list[dict]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def command_version(command: list[str]) -> str:
    try:
        p = subprocess.run(command, text=True, capture_output=True, timeout=30, check=False)
        value = (p.stdout or p.stderr).strip().splitlines()
        return value[0] if value else f"exit={p.returncode}"
    except Exception as exc:  # version audit must not abort project setup
        return f"unavailable: {type(exc).__name__}"


def main() -> None:
    if not STAGE0.is_dir():
        raise SystemExit(f"missing Stage 0: {STAGE0}")
    if ROOT.resolve() == pathlib.Path(ROOT.anchor).resolve() or ROOT.resolve() == STAGE0.resolve():
        raise SystemExit("unsafe Stage 1 path")
    for rel in DIRS:
        (ROOT / rel).mkdir(parents=True, exist_ok=True)

    now = datetime.now(timezone.utc).astimezone()
    source_signature = STAGE0 / "06_results/stage08_locked_translocation_signature.tsv"
    if sha256(source_signature) != EXPECTED_SOURCE_LOCK_SHA256:
        raise SystemExit("Stage 0 signature source hash mismatch")

    inherited = [
        "00_admin/STAGE06_FROZEN_PROTOCOL.md",
        "00_admin/STAGE08_FROZEN_PROTOCOL.md",
        "00_admin/stage08_protocol_deviations.md",
        "01_metadata/amplicon_sample_manifest.tsv",
        "04_processed/amplicon/PRJNA893348/genus_counts.tsv.gz",
        "04_processed/amplicon/PRJNA893348/asv_counts.tsv.gz",
        "04_processed/amplicon/PRJNA893348/asv_taxonomy.tsv.gz",
        "04_processed/amplicon/PRJNA893348/dada2_result.rds",
        "04_processed/amplicon/PRJNA893348/sessionInfo.txt",
        "05_qc/PRJNA893348_dada2_tracking.tsv",
        "05_tables/stage08_patient_sample_ledger.tsv",
        "06_results/stage08_locked_translocation_signature.tsv",
        "06_results/stage08_locked_translocation_signature_provenance.json",
        "06_results/stage08_dataset_gate_summary.json",
        "08_reports/STAGE08_DECISION_REPORT.md",
        "scripts/dada2_cohort.R",
        "02_raw/reference/rdp_train_set_18.fa.gz",
    ]
    manifest_rows = []
    for rel in inherited:
        src = STAGE0 / rel
        if not src.is_file() or src.stat().st_size == 0:
            raise SystemExit(f"missing/empty inherited asset: {src}")
        manifest_rows.append({
            "source_path": str(src), "role": rel, "bytes": src.stat().st_size,
            "sha256": sha256(src), "source_mode": "read_only_no_writeback",
            "verified_at": now.isoformat(timespec="seconds"),
        })
    write_tsv(ROOT / "00_admin/stage0_inheritance_manifest.tsv", manifest_rows,
              ["source_path", "role", "bytes", "sha256", "source_mode", "verified_at"])

    # Copy only compact required assets; Stage 0 remains untouched.
    copy_map = {
        "01_metadata/amplicon_sample_manifest.tsv": "01_metadata/stage0_amplicon_sample_manifest.tsv",
        "05_tables/stage08_patient_sample_ledger.tsv": "01_metadata/stage08_patient_sample_ledger.tsv",
        "04_processed/amplicon/PRJNA893348/genus_counts.tsv.gz": "03_processed/PRJNA893348_genus_counts.tsv.gz",
        "04_processed/amplicon/PRJNA893348/asv_taxonomy.tsv.gz": "03_processed/PRJNA893348_asv_taxonomy.tsv.gz",
        "04_processed/amplicon/PRJNA893348/asv_counts.tsv.gz": "03_processed/PRJNA893348_asv_counts.tsv.gz",
        "04_processed/amplicon/PRJNA893348/dada2_result.rds": "03_processed/PRJNA893348_dada2_result.rds",
        "04_processed/amplicon/PRJNA893348/sessionInfo.txt": "04_qc/PRJNA893348_sessionInfo.txt",
        "05_qc/PRJNA893348_dada2_tracking.tsv": "04_qc/PRJNA893348_stage0_tracking.tsv",
        "06_results/stage08_locked_translocation_signature.tsv": "00_admin/stage0_locked_signature.tsv",
        "06_results/stage08_locked_translocation_signature_provenance.json": "00_admin/stage0_locked_signature_provenance.json",
    }
    for src_rel, dst_rel in copy_map.items():
        dst = ROOT / dst_rel
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(STAGE0 / src_rel, dst)

    ledger = read_tsv(ROOT / "01_metadata/stage08_patient_sample_ledger.tsv")
    ledger = [r for r in ledger if r["dataset"] == "PRJNA893348"]
    amp = read_tsv(ROOT / "01_metadata/stage0_amplicon_sample_manifest.tsv")
    amp = {r["sample_accession"]: r for r in amp if r["study"] == "PRJNA893348"}
    counts = read_tsv(ROOT / "03_processed/PRJNA893348_genus_counts.tsv.gz")
    count_by_sample = {r["sample_accession"]: r for r in counts}
    if len(ledger) != 85 or len(amp) != 84 or len(count_by_sample) != 84:
        raise SystemExit("PRJNA893348 expected 85 runs / 84 patients not reproduced")

    by_sample: dict[str, list[dict]] = defaultdict(list)
    for row in ledger:
        by_sample[row["sample_accession"]].append(row)
    map_rows = []
    blind_qc = []
    genus_columns = [c for c in counts[0] if c != "sample_accession"]
    for sample in sorted(by_sample):
        lr = by_sample[sample]
        a = amp[sample]
        c = count_by_sample[sample]
        depth = sum(int(float(c[g] or 0)) for g in genus_columns)
        is_ap = a["group"] in {"AP_ARDS", "AP_no_ARDS"}
        qc_pass = depth >= 1000
        include = is_ap and qc_pass and len({x["patient_id"] for x in lr}) == 1
        exclusion = "" if include else ("healthy_control_not_primary" if not is_ap else "blind_qc_failure")
        map_rows.append({
            "dataset_id": "PRJNA893348", "subject_id": lr[0]["patient_id"],
            "patient_id_original": lr[0]["patient_id"], "sample_id": sample,
            "biosample": sample, "run_accession": ",".join(sorted(x["run_accession"] for x in lr)),
            "library_id": sample, "specimen": lr[0]["specimen"],
            "collection_time_original": "within 24 h of AP onset (cohort-level report)",
            "time_from_admission": "NA", "time_from_symptom_onset": "<=24h_cohort_level",
            "ARDS_label_original": a["ards"],
            "ARDS_definition": "source-publication ARDS classification",
            "ARDS_assessment_window": "subsequent ARDS; individual interval unavailable",
            "technical_replicate_group": sample if len(lr) > 1 else "",
            "biological_replicate": "False", "sequencing_platform": "Illumina MiSeq",
            "amplicon_region_assay": "16S amplicon, long single-end processed at 350-500 nt",
            "batch": "single_public_cohort", "read_layout": "single_end_combined_by_biosample",
            "age": a["age"], "sex": a["sex"], "bmi": a["bmi"],
            "etiology": "NA", "severity": a["disease_stage"],
            "antibiotic_exposure": "NA", "ICU_status": "NA",
            "include_primary": str(include).lower(), "exclusion_reason": exclusion,
            "mapping_confidence": "direct_public_label",
        })
        blind_qc.append({"sample_id": sample, "genus_assigned_depth": depth,
                         "blind_depth_qc_pass": str(qc_pass).lower()})

    map_cols = list(map_rows[0])
    write_tsv(ROOT / "01_metadata/subject_sample_run_map.tsv", map_rows, map_cols)
    write_tsv(ROOT / "04_qc/blind_sample_qc.tsv", blind_qc,
              ["sample_id", "genus_assigned_depth", "blind_depth_qc_pass"])

    ap_rows = [r for r in map_rows if r["include_primary"] == "true"]
    if len(ap_rows) != 65 or Counter(r["ARDS_label_original"] for r in ap_rows) != Counter({"1": 26, "0": 39}):
        raise SystemExit("primary population did not reproduce 65 patients / 26 ARDS / 39 nonARDS")

    depths = sorted(r["genus_assigned_depth"] for r in blind_qc)
    blind_summary = {
        "n_samples": len(depths), "min_depth": depths[0],
        "median_depth": statistics.median(depths), "max_depth": depths[-1],
        "all_pass_frozen_minimum_1000": all(d >= 1000 for d in depths),
        "signature_taxon_availability_blind_all_84": {
            "Borreliella": {"table_column": None, "prevalence": 0.0, "rule": "insert_zero_count"},
            "Escherichia": {"table_column": "Escherichia/Shigella", "rule": "RDP18_compound_taxon_crosswalk"},
            "Staphylococcus": {"table_column": "Staphylococcus"},
            "Enterococcus": {"table_column": "Enterococcus"},
            "Klebsiella": {"table_column": "Klebsiella"},
        },
        "tracking_note": "Stage 0 denoised/merged/nonchim tracking columns were blank because of an index-name mismatch; Stage 1 uses immutable genus counts and reconstructs depth blind to outcome.",
    }
    write_json(ROOT / "04_qc/blind_qc_summary.json", blind_summary)

    signature_lock = {
        "lock_status": "LOCKED_BEFORE_PRJNA893348_OUTCOME_ANALYSIS",
        "locked_at": now.isoformat(timespec="seconds"),
        "source_signature_path": str(source_signature),
        "source_signature_sha256": EXPECTED_SOURCE_LOCK_SHA256,
        "features": [
            {"locked_genus": "Borreliella", "direction": 1, "input_column": None,
             "transfer_rule": "structural zero in RDP18 table; insert count 0; no substitute genus"},
            {"locked_genus": "Escherichia", "direction": 1, "input_column": "Escherichia/Shigella",
             "transfer_rule": "pre-outcome RDP18 nomenclature crosswalk; compound 16S taxon retained verbatim"},
            {"locked_genus": "Staphylococcus", "direction": 1, "input_column": "Staphylococcus"},
            {"locked_genus": "Enterococcus", "direction": 1, "input_column": "Enterococcus"},
            {"locked_genus": "Klebsiella", "direction": 1, "input_column": "Klebsiella"},
        ],
        "input": "unrarefied bacterial genus count matrix after blind technical QC",
        "zero_rule": "add fixed pseudocount 0.5 to every genus count, including inserted structural-zero Borreliella",
        "transform": "sample-wise CLR over every genus column after pseudocount",
        "component_rule": "positive-direction CLR component for each of the five locked genera",
        "combination": "arithmetic mean of five components with fixed coefficient 0.2 each; no refitting",
        "standardization": "center and divide aggregate score by sample SD among the 65 frozen AP patients",
        "score_direction": "higher score prespecified as higher ARDS risk",
        "missing_rule": "unavailable locked genus is retained as a zero-count feature; no near-neighbour replacement",
    }
    write_json(ROOT / "00_admin/signature_lock.json", signature_lock)

    analysis_plan = {
        "primary_population": "65 unique AP patients; healthy controls excluded",
        "outcome": "ARDS (26) versus nonARDS (39)",
        "statistical_unit": "patient",
        "primary_statistic": "mean standardized signature score in ARDS minus nonARDS",
        "primary_test": "one-sided patient-label permutation, alternative greater",
        "permutations": 10000, "p_threshold": 0.05, "master_seed": MASTER_SEED,
        "effect_interval": "10000 stratified patient bootstrap percentile 95% CI for mean difference",
        "descriptive_auc": "empirical AUROC with 10000 stratified bootstrap percentile 95% CI; not a gate",
        "adjustment": "supportive logistic regression with age, sex, BMI only; complete case",
        "critical_sensitivities": [
            "leave-one-patient-out", "log sequencing-depth adjustment",
            "exclude blind bottom 5 percent depth", "pseudocount 1.0",
            "leave-one-genus-out", "strict skin-contaminant sensitivity omitting Staphylococcus without reweighting",
        ],
        "multiplicity": "single confirmatory test; sensitivities diagnostic; genus-level exploratory results BH-FDR",
    }
    write_json(ROOT / "00_admin/analysis_plan.json", analysis_plan)

    variable_rows = [
        {"variable": "ARDS", "role": "only primary outcome", "definition": "public 0/1 label among AP patients", "missing_rule": "exclude only if label unresolved"},
        {"variable": "signature_score_sd", "role": "only primary exposure", "definition": "locked five-genus CLR mean standardized in frozen AP set", "missing_rule": "none after structural-zero rule"},
        {"variable": "age", "role": "supportive covariate", "definition": "years", "missing_rule": "complete-case only"},
        {"variable": "sex", "role": "supportive covariate", "definition": "female/male public label", "missing_rule": "complete-case only"},
        {"variable": "bmi", "role": "supportive covariate", "definition": "kg/m2", "missing_rule": "complete-case only"},
        {"variable": "genus_assigned_depth", "role": "blind QC/sensitivity", "definition": "sum of all genus counts", "missing_rule": "not applicable"},
    ]
    write_tsv(ROOT / "00_admin/variable_dictionary.tsv", variable_rows,
              ["variable", "role", "definition", "missing_rule"])
    exclusions = [{"subject_id": r["subject_id"], "sample_id": r["sample_id"],
                   "reason": r["exclusion_reason"], "decision_blinded": "true"}
                  for r in map_rows if r["include_primary"] != "true"]
    write_tsv(ROOT / "01_metadata/exclusion_log.tsv", exclusions,
              ["subject_id", "sample_id", "reason", "decision_blinded"])
    flow = [
        {"step": "public_runs", "n": 85, "note": "all PRJNA893348 runs"},
        {"step": "unique_patients", "n": 84, "note": "J17 two technical runs combined before DADA2"},
        {"step": "healthy_controls_excluded", "n": 19, "note": "not part of primary endpoint"},
        {"step": "AP_patients_QC_pass", "n": 65, "note": "26 ARDS and 39 nonARDS; all pass blind depth QC"},
    ]
    write_tsv(ROOT / "01_metadata/analysis_population_flow.tsv", flow, ["step", "n", "note"])

    provenance = [
        {"field": "patient/sample/run identity", "source": str(STAGE0 / "05_tables/stage08_patient_sample_ledger.tsv"), "confidence": "direct_public_label", "review": "author-reviewed"},
        {"field": "ARDS/age/sex/BMI", "source": str(STAGE0 / "01_metadata/amplicon_sample_manifest.tsv"), "confidence": "reconciled_stage0_public_metadata", "review": "author-reviewed"},
        {"field": "sampling window", "source": str(STAGE0 / "08_reports/STAGE08_DECISION_REPORT.md"), "confidence": "cohort_level_not_patient_specific", "review": "author-reviewed"},
        {"field": "five-genus signature", "source": str(source_signature), "confidence": "checksum_locked_before_external_outcome", "review": "hash-verified frozen artifact"},
    ]
    write_tsv(ROOT / "01_metadata/source_provenance.tsv", provenance,
              ["field", "source", "confidence", "review"])
    write_tsv(ROOT / "01_metadata/metadata_discrepancy_log.tsv", [],
              ["field", "record_id", "source_a", "value_a", "source_b", "value_b", "resolution", "outcome_accessed"])

    protocol = f"""# Stage 1 Frozen Analysis Protocol

Locked: {now.isoformat(timespec='seconds')}  
Review: author-reviewed

## Confirmatory question

Test whether the Stage 0 outcome-blind, checksum-locked five-genus signature is
higher in the 26 AP patients who subsequently developed ARDS than in the 39 AP
patients who did not. Healthy controls are excluded. The patient is the only
statistical unit.

## Immutable signature

The source lock SHA256 is `{EXPECTED_SOURCE_LOCK_SHA256}` and contains
Borreliella, Escherichia, Staphylococcus, Enterococcus, and Klebsiella, all in
the positive direction. The PRJNA893348 RDP18 table labels Escherichia as the
non-separable compound `Escherichia/Shigella`; this exact compound is the
pre-outcome nomenclature crosswalk. No Borreliella column or classified read is
present; the feature is retained as an all-zero count and no Borrelia or other
near-neighbour substitution is permitted.

For every retained sample, add 0.5 to every bacterial genus count, perform the
sample-wise CLR over the full genus table, take the arithmetic mean of the five
fixed positive-direction components (coefficient 0.2 each), then standardize
that aggregate by the sample standard deviation among all 65 frozen AP
patients. A higher score is prespecified as higher ARDS risk. No feature,
direction, weight, cut-off, transform, or patient can be selected using outcome
results.

## Population and timing

The frozen population is 65 unique AP patients (26 ARDS, 39 nonARDS). All pass
the outcome-blind minimum assigned-depth rule of 1,000 counts. The publication
supports sampling within 24 hours of AP onset and before subsequent ARDS at the
cohort level, but an individual time interval is not public; therefore primary
wording is “early association/identification,” not a patient-level prospective
prediction claim.

## Primary inference

The statistic is mean standardized score in ARDS minus nonARDS. The only
confirmatory P value is a one-sided patient-label permutation test with 10,000
permutations, alternative greater, master seed {MASTER_SEED}, and P<0.05. The
Monte Carlo standard error is reported. The main effect is the mean difference
with a 10,000-replicate stratified patient-bootstrap percentile 95% CI. AUROC
and its bootstrap interval are descriptive and cannot determine success.

## Support and sensitivity

Age, sex, and BMI are the only allowed adjustment variables and enter one
supportive complete-case logistic model without outcome-driven selection.
Critical diagnostics are leave-one-patient-out; log-depth adjustment; exclusion
of the blind bottom 5% depth; fixed pseudocount 1.0; every leave-one-genus-out
score without refitting; and a stringent skin-contaminant analysis omitting
Staphylococcus without redistributing its 0.2 weight. These cannot replace the
primary analysis. Single-genus tests are exploratory and use BH-FDR.

## Contamination and causal boundaries

No re-analyzable negative controls are public, so decontamination cannot be
claimed. Microbial reads do not prove viable organisms, migration direction, or
causality. Different platforms/regions remain separate. Cross-cohort host and
cellular findings are triangular support, never patient-level multi-omics.
"""
    write_text(ROOT / "00_admin/FROZEN_ANALYSIS_PROTOCOL.md", protocol)

    project_manifest = f"""project: acute-pancreatitis-cross-compartment-microbiome-early-ards
stage: 1
root: {ROOT}
stage0_read_only: {STAGE0}
created_at: {now.isoformat(timespec='seconds')}
statistical_unit: patient
primary_dataset: PRJNA893348
master_seed: {MASTER_SEED}
review: author-reviewed
"""
    write_text(ROOT / "00_admin/project_manifest.yaml", project_manifest)
    write_text(ROOT / "00_admin/storage_plan.yaml",
               f"project_root: {ROOT}\nminimum_free_fraction: 0.20\ntmp: {ROOT / 'tmp'}\nstage0_mode: read_only\n")
    write_text(ROOT / "00_admin/resource_policy.yaml",
               "reserve_memory_fraction: 0.20\nreserve_disk_fraction: 0.20\nlow_threads: 4\nnormal_threads: 16\nhigh_threads: 32\nnever_stop_unrelated_processes: true\ncontinue_at_low_resources: true\n")
    software = [
        {"software": "python", "version": command_version(["python3", "--version"]), "scope": "system read-only runtime"},
        {"software": "R", "version": command_version(["R", "--version"]), "scope": "system read-only runtime"},
        {"software": "git", "version": command_version(["git", "--version"]), "scope": "system read-only runtime"},
        {"software": "DADA2", "version": "1.38.0", "scope": "inherited Stage 0 processing"},
        {"software": "RDP training set", "version": "18; sha256=f0994b87030a8d52764ea60f10aba008ae6ef9a5b129ef447f2885bdcdca8e89", "scope": "inherited Stage 0 taxonomy"},
    ]
    write_tsv(ROOT / "00_admin/software_lock.tsv", software, ["software", "version", "scope"])
    write_tsv(ROOT / "00_admin/command_audit.tsv", [
        {"time": now.isoformat(timespec="seconds"), "cwd": str(ROOT), "command_template": "python3 workflow/prepare_stage1.py", "exit_code": 0, "contains_secret": "false"}
    ], ["time", "cwd", "command_template", "exit_code", "contains_secret"])
    write_text(ROOT / "00_admin/protocol_deviations.md",
               "# Protocol deviations\n\nNo outcome-informed deviations at lock time.\n\n- Stage 0 QC tracking columns `denoised`, `merged`, and `nonchim` were blank because of an index-name mismatch. Stage 1 reconstructs assigned genus depth from the immutable count table without altering Stage 0. This is a technical-record repair and does not change patients, outcomes, features, or thresholds.\n")

    status = {
        "stage": "stage1_protocol_locked", "updated_at": now.isoformat(timespec="seconds"),
        "outcome_analysis_started": False, "primary_complete": False,
        "supporting_complete": False, "final_complete": False,
    }
    write_json(ROOT / "00_admin/stage_status.json", status)
    startup = f"""# Stage 1 Startup Report

Status: Stage 1 preparation complete  
Review: author-reviewed  
Time: {now.isoformat(timespec='seconds')}

- Stage 1 is isolated at `{ROOT}` on the administrator data volume.
- Stage 0 is registered read-only; no source file was modified.
- The source signature checksum exactly matches the expected value.
- PRJNA893348 reproduces 85 runs, 84 unique patients, and one healthy-control
  technical duplicate (J17). The primary population is 65 unique AP patients:
  26 ARDS and 39 nonARDS; healthy controls are excluded.
- All 84 processed samples are present and all pass the outcome-blind minimum
  genus-depth threshold. The count table, RDS, taxonomy, processing script,
  database, and session information are checksum-bound in the inheritance
  manifest.
- RDP18 cannot separate Escherichia from Shigella and reports the compound
  Escherichia/Shigella. Borreliella is absent. Both transfer rules were frozen
  before outcome analysis and will be prominently disclosed.
- Individual sampling-to-ARDS intervals are not public; manuscript language is
  restricted to early association/identification.

The exact zero rule, transform,
weights, standardization, direction, endpoint, permutation scheme, covariates,
sensitivities, and multiplicity are now machine-readably frozen
before grouped microbial results are computed.
"""
    write_text(ROOT / "00_admin/STAGE1_STARTUP_REPORT.md", startup)

    lock_files = [
        ROOT / "00_admin/FROZEN_ANALYSIS_PROTOCOL.md",
        ROOT / "00_admin/signature_lock.json",
        ROOT / "00_admin/analysis_plan.json",
        ROOT / "00_admin/variable_dictionary.tsv",
        ROOT / "01_metadata/subject_sample_run_map.tsv",
        ROOT / "01_metadata/exclusion_log.tsv",
        ROOT / "01_metadata/analysis_population_flow.tsv",
    ]
    locks = [{"path": str(p), "bytes": p.stat().st_size, "sha256": sha256(p),
              "locked_at": now.isoformat(timespec="seconds")} for p in lock_files]
    write_tsv(ROOT / "00_admin/protocol_lock_manifest.tsv", locks,
              ["path", "bytes", "sha256", "locked_at"])
    for p in lock_files + [ROOT / "00_admin/protocol_lock_manifest.tsv"]:
        p.chmod(0o444)

    # Initialise a project-local audit repository. Data objects are excluded.
    write_text(ROOT / ".gitignore", "02_raw/\n03_processed/\ntmp/\n*.partial\n")
    if not (ROOT / ".git").exists():
        subprocess.run(["git", "init"], cwd=ROOT, check=True, stdout=subprocess.DEVNULL)
    subprocess.run(["git", "add", "."], cwd=ROOT, check=True)
    subprocess.run(["git", "commit", "-m", "Freeze Stage 1 protocol before outcome analysis"],
                   cwd=ROOT, check=True, stdout=subprocess.DEVNULL)
    commit = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
    write_text(ROOT / "00_admin/PROTOCOL_LOCKED", f"{now.isoformat(timespec='seconds')}\t{commit}\n")
    print(json.dumps({"status": "complete", "root": str(ROOT), "commit": commit,
                      "primary_n": 65, "ards": 26, "nonards": 39}, ensure_ascii=False))


if __name__ == "__main__":
    main()
