#!/usr/bin/env python3
"""Resumable ENA download and checksum validation for PRJNA428535."""

import csv
import hashlib
import json
import os
import pathlib
import subprocess
from datetime import datetime, timezone

ROOT = pathlib.Path(os.environ.get("STAGE1_ROOT", pathlib.Path.cwd()))
MANIFEST = ROOT / "01_metadata/PRJNA428535_ena_download_manifest.tsv"
OUT = ROOT / "02_raw/PRJNA428535"


def md5(path):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(4 * 1024 * 1024), b""):
            h.update(b)
    return h.hexdigest()


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    rows = list(csv.DictReader(open(MANIFEST, encoding="utf-8"), delimiter="\t"))
    if len(rows) != 124:
        raise SystemExit("ENA manifest does not contain 124 runs")
    aria = ROOT / "tmp/PRJNA428535_aria2.txt"
    with open(aria, "w", encoding="utf-8", newline="\n") as fh:
        for r in rows:
            urls = r["fastq_ftp"].split(";")
            if len(urls) != 1:
                raise SystemExit(f"unexpected multifile run {r['run_accession']}")
            name = pathlib.PurePosixPath(urls[0]).name
            fh.write("https://" + urls[0] + "\n")
            fh.write(f"  dir={OUT}\n  out={name}\n")
    cmd = ["aria2c", "--input-file", str(aria), "--continue=true",
           "--max-concurrent-downloads=4", "--split=4", "--max-connection-per-server=4",
           "--min-split-size=5M", "--file-allocation=none", "--auto-file-renaming=false",
           "--allow-overwrite=false", "--summary-interval=60", "--console-log-level=notice"]
    subprocess.run(cmd, check=True)
    checks = []
    failures = []
    for r in rows:
        url = r["fastq_ftp"].split(";")[0]
        p = OUT / pathlib.PurePosixPath(url).name
        got = md5(p) if p.is_file() else "MISSING"
        ok = got == r["fastq_md5"]
        checks.append({"run_accession": r["run_accession"], "path": str(p),
                       "expected_md5": r["fastq_md5"], "observed_md5": got,
                       "bytes": p.stat().st_size if p.is_file() else 0, "verified": ok})
        if not ok:
            failures.append(r["run_accession"])
    out = ROOT / "04_qc/PRJNA428535_download_checksums.tsv"
    with open(out.with_suffix(".tsv.partial"), "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(checks[0]), delimiter="\t")
        w.writeheader(); w.writerows(checks)
    os.replace(out.with_suffix(".tsv.partial"), out)
    if failures:
        (ROOT / "logs/PRJNA428535_DOWNLOAD_FAILED").write_text("\n".join(failures) + "\n")
        raise SystemExit(f"checksum failures: {len(failures)}")
    (ROOT / "logs/PRJNA428535_DOWNLOAD_COMPLETE").write_text(
        datetime.now(timezone.utc).astimezone().isoformat(timespec="seconds") + "\n")
    print(json.dumps({"runs": len(rows), "verified": len(checks), "status": "COMPLETE"}))


if __name__ == "__main__":
    main()
