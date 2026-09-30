from __future__ import annotations

import csv
import gzip
import json
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def read_tsv(path: Path) -> list[dict[str, str]]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


class RepositoryTests(unittest.TestCase):
    def test_required_files(self) -> None:
        required = [
            "README.md", "LICENSE", "CITATION.cff", "requirements.txt", "environment.yml",
            "scripts/analysis/reproduce_primary.py", "scripts/analysis/reproduce_secondary.py",
            "scripts/analysis/verify_reproduction.py", "scripts/figures/make_figures.R",
            "data/processed/PRJNA893348_genus_counts.tsv.gz",
            "data/processed/PRJNA428535_genus_counts.tsv.gz",
            "data/metadata/subject_sample_run_map.tsv",
            "data/metadata/stage08_patient_sample_ledger.tsv",
        ]
        for relative in required:
            path = ROOT / relative
            self.assertTrue(path.is_file() and path.stat().st_size > 0, relative)

    def test_primary_population(self) -> None:
        rows = [r for r in read_tsv(ROOT / "data/metadata/subject_sample_run_map.tsv") if r["include_primary"] == "true"]
        self.assertEqual(len(rows), 65)
        self.assertEqual(sum(int(r["ARDS_label_original"]) for r in rows), 26)
        self.assertEqual(sum(1 - int(r["ARDS_label_original"]) for r in rows), 39)
        self.assertEqual(len({r["subject_id"] for r in rows}), 65)

    def test_paired_and_host_counts(self) -> None:
        ledger = [r for r in read_tsv(ROOT / "data/metadata/stage08_patient_sample_ledger.tsv") if r["dataset"] == "PRJNA428535"]
        self.assertEqual(len(ledger), 124)
        self.assertEqual(len({r["patient_id"] for r in ledger}), 62)
        host = read_tsv(ROOT / "data/host_response/host_module_scores.tsv")
        self.assertEqual(len(host), 87)

    def test_primary_result_values(self) -> None:
        result = json.loads((ROOT / "results/analysis/primary_test.json").read_text(encoding="utf-8"))
        self.assertAlmostEqual(result["observed_mean_difference_SD"], 0.5947436540143409, places=12)
        self.assertAlmostEqual(result["bootstrap_95CI"][0], 0.11441789386335492, places=12)
        self.assertAlmostEqual(result["bootstrap_95CI"][1], 1.0653136209549643, places=12)
        self.assertAlmostEqual(result["one_sided_permutation_p"], 0.009299070092990702, places=15)

    def test_no_machine_paths_or_credentials(self) -> None:
        banned = [
            "/data/" + "admin/",
            "c:\\" + "users\\" + "aaa12",
            "master" + "2333",
            "begin openssh "
            + "private key",
            "begin rsa "
            + "private key",
        ]
        suffixes = {".py", ".r", ".sh", ".md", ".json", ".tsv", ".txt", ".yml", ".yaml", ".cff"}
        for path in ROOT.rglob("*"):
            if not path.is_file() or path.suffix.lower() not in suffixes:
                continue
            text = path.read_text(encoding="utf-8", errors="ignore").lower()
            for token in banned:
                self.assertNotIn(token, text, f"{token} in {path.relative_to(ROOT)}")


if __name__ == "__main__":
    unittest.main()
