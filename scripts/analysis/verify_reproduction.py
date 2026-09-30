#!/usr/bin/env python3
"""Compare regenerated outputs with frozen publication-candidate results."""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
GENERATED = ROOT / "results" / "analysis"
PUBLISHED = ROOT / "data" / "published_results"
HOST = ROOT / "data" / "host_response"
ATOL = 1e-8
RTOL = 1e-9


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def as_float(value: str) -> float | None:
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def compare_tsv(actual: Path, expected: Path, atol: float = ATOL, rtol: float = RTOL) -> list[str]:
    problems: list[str] = []
    if not actual.is_file() or not expected.is_file():
        return [f"missing comparison file: {actual if not actual.is_file() else expected}"]
    a_rows, e_rows = read_tsv(actual), read_tsv(expected)
    if len(a_rows) != len(e_rows):
        return [f"{actual.name}: row count {len(a_rows)} != {len(e_rows)}"]
    if a_rows and e_rows and list(a_rows[0]) != list(e_rows[0]):
        problems.append(f"{actual.name}: columns differ")
        return problems
    for index, (a_row, e_row) in enumerate(zip(a_rows, e_rows), start=2):
        for column in e_row:
            a_value, e_value = a_row[column], e_row[column]
            a_number, e_number = as_float(a_value), as_float(e_value)
            if a_number is not None and e_number is not None:
                if not math.isclose(a_number, e_number, rel_tol=rtol, abs_tol=atol):
                    problems.append(f"{actual.name}:{index}:{column}: {a_number} != {e_number}")
            elif a_value != e_value:
                problems.append(f"{actual.name}:{index}:{column}: {a_value!r} != {e_value!r}")
            if len(problems) >= 25:
                return problems
    return problems


def main() -> int:
    comparisons = [
        "primary_effects.tsv",
        "patient_signature_components.tsv",
        "leave_one_patient_out.tsv",
        "contamination_sensitivity.tsv",
        "covariate_support.tsv",
        "primary_component_details.tsv",
        "sensitivity_with_ci.tsv",
        "alternative_compositional_sensitivity.tsv",
        "clinical_baseline_by_ards.tsv",
        "clinical_missingness.tsv",
        "firth_logistic_models.tsv",
        "patient_outcome_timing.tsv",
        "paired_scores.tsv",
        "paired_results.tsv",
        "genus_pair_concordance.tsv",
        "paired_correlation_ci.tsv",
        "paired_genus_effects.tsv",
        "host_module_dictionary.tsv",
        "host_module_overlap.tsv",
    ]
    problems: list[str] = []
    for name in comparisons:
        # Firth estimates can vary at the last optimizer digits across SciPy
        # builds; primary effect and permutation results retain the strict gate.
        if name == "firth_logistic_models.tsv":
            problems.extend(compare_tsv(GENERATED / name, PUBLISHED / name, atol=1e-6, rtol=1e-7))
        else:
            problems.extend(compare_tsv(GENERATED / name, PUBLISHED / name))
    for name in ["host_module_results.tsv", "host_module_scores.tsv"]:
        problems.extend(compare_tsv(HOST / name, PUBLISHED / name))

    primary = json.loads((GENERATED / "primary_test.json").read_text(encoding="utf-8"))
    expected = {
        "n": 65,
        "ARDS": 26,
        "nonARDS": 39,
        "observed_mean_difference_SD": 0.5947436540143409,
        "one_sided_permutation_p": 0.009299070092990702,
    }
    for key, value in expected.items():
        actual = primary.get(key)
        if isinstance(value, float):
            if actual is None or not math.isclose(float(actual), value, rel_tol=RTOL, abs_tol=ATOL):
                problems.append(f"primary_test.json:{key}: {actual} != {value}")
        elif actual != value:
            problems.append(f"primary_test.json:{key}: {actual} != {value}")

    sensitivity = read_tsv(GENERATED / "sensitivity_with_ci.tsv")[0]
    sensitivity_p = float(sensitivity["one_sided_permutation_p"])
    if not math.isclose(sensitivity_p, 0.008799120087991202, rel_tol=0.0, abs_tol=1e-15):
        problems.append("secondary permutation-seed sensitivity P value changed")
    if math.isclose(sensitivity_p, float(primary["one_sided_permutation_p"]), rel_tol=0.0, abs_tol=1e-15):
        problems.append("primary and sensitivity P values were incorrectly collapsed")

    required_figures = list((ROOT / "results" / "figures").glob("Figure_[1-6]_*.pdf"))
    if len(required_figures) != 6 or any(path.stat().st_size == 0 for path in required_figures):
        problems.append("six non-empty final PDF figures are required")
    source_tables = list((ROOT / "data" / "figure_source_data").glob("*.tsv"))
    if len(source_tables) != 16:
        problems.append(f"expected 16 figure source-data tables, found {len(source_tables)}")

    report = {
        "status": "PASS" if not problems else "FAIL",
        "numeric_tolerance": {
            "default_absolute": ATOL,
            "default_relative": RTOL,
            "firth_optimizer_absolute": 1e-6,
            "firth_optimizer_relative": 1e-7,
        },
        "compared_result_tables": comparisons,
        "host_tables_checked": ["host_module_results.tsv", "host_module_scores.tsv"],
        "primary_result": expected,
        "secondary_seed_sensitivity_p": sensitivity_p,
        "figure_pdf_count": len(required_figures),
        "figure_source_table_count": len(source_tables),
        "problems": problems,
    }
    output = GENERATED / "reproduction_report.json"
    output.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2, ensure_ascii=False))
    return 0 if not problems else 1


if __name__ == "__main__":
    raise SystemExit(main())
