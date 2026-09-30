#!/usr/bin/env python3
"""Paired blood-neutrophil support analysis for PRJNA428535."""

import csv
import gzip
import json
import math
import os
import pathlib
import shutil
from collections import defaultdict
from datetime import datetime, timezone

import numpy as np
from scipy import stats

ROOT = pathlib.Path(os.environ.get("STAGE1_ROOT", pathlib.Path.cwd()))
STAGE0 = pathlib.Path(os.environ.get("STAGE0_ROOT", ROOT.parent / "gut_pancreas_twin_stage0"))
SEED = 20260929
B = 10_000
FEATURES = ["Borreliella", "Escherichia", "Staphylococcus", "Enterococcus", "Klebsiella"]
INPUT_MAP = {"Borreliella": None, "Escherichia": "Escherichia/Shigella",
             "Staphylococcus": "Staphylococcus", "Enterococcus": "Enterococcus", "Klebsiella": "Klebsiella"}


def read_tsv(path):
    path = pathlib.Path(path); op = gzip.open if path.suffix == ".gz" else open
    with op(path, "rt", encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write_text(path, text):
    path = pathlib.Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".partial"); tmp.write_text(text, encoding="utf-8", newline="\n"); os.replace(tmp, path)


def write_json(path, obj):
    write_text(path, json.dumps(obj, ensure_ascii=False, indent=2) + "\n")


def write_tsv(path, rows, cols=None):
    path = pathlib.Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    cols = cols or (list(rows[0]) if rows else [])
    tmp = path.with_suffix(path.suffix + ".partial")
    with tmp.open("w", encoding="utf-8", newline="") as fh:
        w=csv.DictWriter(fh,fieldnames=cols,delimiter="\t",extrasaction="ignore"); w.writeheader(); w.writerows(rows)
    os.replace(tmp,path)


def make_scores(counts, genera, omitted=None):
    omitted = omitted or set(); g=list(genera); x=counts.astype(float)
    if "Borreliella" not in g:
        g.append("Borreliella"); x=np.column_stack([x,np.zeros(x.shape[0])])
    log=np.log(x+0.5); clr=log-log.mean(axis=1,keepdims=True)
    comp={}; total=np.zeros(x.shape[0])
    for f in FEATURES:
        col=INPUT_MAP[f] or "Borreliella"; comp[f]=clr[:,g.index(col)]
        if f not in omitted: total += 0.2*comp[f]
    total=(total-total.mean())/total.std(ddof=1)
    return total,comp


def corr_perm(a,b,seed):
    rng=np.random.default_rng(seed); obs=stats.spearmanr(a,b).statistic
    vals=np.asarray([stats.spearmanr(a,rng.permutation(b)).statistic for _ in range(B)])
    return float(obs),float((1+np.sum(vals>=obs))/(B+1))


def paired_diff(a,b,seed):
    # b minus a, interpreted as neutrophil minus blood.
    d=b-a; obs=float(d.mean()); rng=np.random.default_rng(seed)
    perm=np.asarray([np.mean(d*rng.choice([-1,1],len(d))) for _ in range(B)])
    p=float((1+np.sum(np.abs(perm)>=abs(obs)))/(B+1))
    boot=np.asarray([np.mean(rng.choice(d,len(d),True)) for _ in range(B)])
    return obs,np.quantile(boot,[.025,.975]),p


def main():
    outdir=ROOT/"06_supporting/blood_neutrophil"; outdir.mkdir(parents=True,exist_ok=True)
    ledger_src=STAGE0/"05_tables/stage08_patient_sample_ledger.tsv"
    ledger_dst=ROOT/"01_metadata/stage08_patient_sample_ledger.tsv"
    if not ledger_dst.exists(): shutil.copy2(ledger_src,ledger_dst)
    ledger=[r for r in read_tsv(ledger_dst) if r["dataset"]=="PRJNA428535"]
    counts_rows=read_tsv(ROOT/"03_processed/PRJNA428535/genus_counts.tsv.gz")
    genera=[x for x in counts_rows[0] if x!="run_accession"]
    cb={r["run_accession"]:r for r in counts_rows}
    if len(ledger)!=124 or set(cb)!=set(r["run_accession"] for r in ledger): raise RuntimeError("run mapping mismatch")
    ledger.sort(key=lambda r:r["run_accession"])
    runs=[r["run_accession"] for r in ledger]
    x=np.asarray([[int(float(cb[run][g] or 0)) for g in genera] for run in runs],float)
    depth=x.sum(axis=1)
    score,comp=make_scores(x,genera)
    strict,_=make_scores(x,genera,omitted={"Staphylococcus"})
    by_patient=defaultdict(dict)
    for i,r in enumerate(ledger):
        by_patient[r["patient_id"]][r["specimen"]]=(i,r)
    if len(by_patient)!=62 or any(set(v)!={"blood","neutrophils"} for v in by_patient.values()): raise RuntimeError("pairing failure")

    pair_rows=[]
    for pid in sorted(by_patient):
        ib,rb=by_patient[pid]["blood"]; in_,rn=by_patient[pid]["neutrophils"]
        row={"patient_id":pid,"group":rb["group"],"blood_run":runs[ib],"neutrophil_run":runs[in_],
             "blood_depth":int(depth[ib]),"neutrophil_depth":int(depth[in_]),
             "blood_score":float(score[ib]),"neutrophil_score":float(score[in_]),
             "blood_strict_score":float(strict[ib]),"neutrophil_strict_score":float(strict[in_]),
             "pair_qc_ge100":str(depth[ib]>=100 and depth[in_]>=100).lower(),
             "pair_qc_ge1000":str(depth[ib]>=1000 and depth[in_]>=1000).lower()}
        for f in FEATURES:
            row[f"blood_{f}_detected"]=int((x[ib,genera.index(INPUT_MAP[f])] if INPUT_MAP[f] in genera else 0)>0)
            row[f"neutrophil_{f}_detected"]=int((x[in_,genera.index(INPUT_MAP[f])] if INPUT_MAP[f] in genera else 0)>0)
        pair_rows.append(row)
    write_tsv(outdir/"paired_scores.tsv",pair_rows)
    keep=np.asarray([r["pair_qc_ge100"]=="true" for r in pair_rows])
    a=np.asarray([r["blood_score"] for r in pair_rows])[keep]
    b=np.asarray([r["neutrophil_score"] for r in pair_rows])[keep]
    rho,p_corr=corr_perm(a,b,SEED+901)
    d,ci,p_diff=paired_diff(a,b,SEED+902)
    results=[{"analysis":"all_pairs_ge100","n_pairs":int(keep.sum()),"spearman_rho":rho,"pairing_permutation_p":p_corr,
              "neutrophil_minus_blood_mean":d,"paired_bootstrap_CI_low":float(ci[0]),"paired_bootstrap_CI_high":float(ci[1]),"signflip_p_two_sided":p_diff}]
    keep1000=np.asarray([r["pair_qc_ge1000"]=="true" for r in pair_rows])
    if keep1000.sum()>=10:
        aa=np.asarray([r["blood_score"] for r in pair_rows])[keep1000]; bb=np.asarray([r["neutrophil_score"] for r in pair_rows])[keep1000]
        rr,pp=corr_perm(aa,bb,SEED+903); dd,cc,pd=paired_diff(aa,bb,SEED+904)
        results.append({"analysis":"pairs_ge1000_sensitivity","n_pairs":int(keep1000.sum()),"spearman_rho":rr,"pairing_permutation_p":pp,
                        "neutrophil_minus_blood_mean":dd,"paired_bootstrap_CI_low":float(cc[0]),"paired_bootstrap_CI_high":float(cc[1]),"signflip_p_two_sided":pd})
    aa=np.asarray([r["blood_strict_score"] for r in pair_rows])[keep]; bb=np.asarray([r["neutrophil_strict_score"] for r in pair_rows])[keep]
    rr,pp=corr_perm(aa,bb,SEED+905); dd,cc,pd=paired_diff(aa,bb,SEED+906)
    results.append({"analysis":"omit_Staphylococcus_fixed_weight","n_pairs":int(keep.sum()),"spearman_rho":rr,"pairing_permutation_p":pp,
                    "neutrophil_minus_blood_mean":dd,"paired_bootstrap_CI_low":float(cc[0]),"paired_bootstrap_CI_high":float(cc[1]),"signflip_p_two_sided":pd})
    write_tsv(outdir/"paired_results.tsv",results)

    concord=[]
    for f in FEATURES:
        bv=np.asarray([r[f"blood_{f}_detected"] for r in pair_rows],int); nv=np.asarray([r[f"neutrophil_{f}_detected"] for r in pair_rows],int)
        concord.append({"genus":f,"n_pairs":len(pair_rows),"blood_detected":int(bv.sum()),"neutrophil_detected":int(nv.sum()),
                        "both_detected":int(np.sum((bv==1)&(nv==1))),"both_absent":int(np.sum((bv==0)&(nv==0))),
                        "agreement_fraction":float(np.mean(bv==nv))})
    write_tsv(outdir/"genus_pair_concordance.tsv",concord)

    subgroup=[]
    for group in sorted(set(r["group"] for r in pair_rows)):
        idx=np.asarray([r["group"]==group and r["pair_qc_ge100"]=="true" for r in pair_rows])
        av=np.asarray([r["blood_score"] for r in pair_rows])[idx]; bv=np.asarray([r["neutrophil_score"] for r in pair_rows])[idx]
        subgroup.append({"group":group,"n_pairs":int(idx.sum()),"spearman_rho":float(stats.spearmanr(av,bv).statistic),
                         "mean_neutrophil_minus_blood":float(np.mean(bv-av))})
    write_tsv(outdir/"paired_group_descriptives.tsv",subgroup)
    summary={"dataset":"PRJNA428535","runs":124,"patients":62,"all_pairs_primary_support":results[0],
             "sensitivities":results[1:],"genus_concordance":concord,
             "limitations":["low-biomass 16S","no public sequenced blanks","Ion Torrent V3 amplicon","cannot establish viability or intracellular carriage or migration direction"],
             "role":"within-person cellular-compartment support only","review":"author-reviewed"}
    write_json(outdir/"blood_neutrophil_summary.json",summary)
    report=f"""# Paired blood–neutrophil support

Review: author-reviewed

All 124 runs were checksum-verified and mapped to 62 unique subject pairs. After
the frozen technical threshold, {int(keep.sum())} pairs entered the primary support
analysis. The patient-level blood–neutrophil score correlation was rho={rho:.3f}
(pairing-permutation P={p_corr:.4g}). The mean neutrophil-minus-blood difference
was {d:.3f} SD (paired bootstrap 95%CI {ci[0]:.3f} to {ci[1]:.3f}; two-sided
sign-flip P={p_diff:.4g}).

This module tests compartment concordance and enrichment only. It is low-biomass
Ion Torrent V3 16S data without public sequenced blanks; qPCR controls reported by
the source paper do not prove viable bacteria, intracellular carriage, or migration
direction. It is not an ARDS validation cohort.
"""
    write_text(outdir/"BLOOD_NEUTROPHIL_REPORT.md",report)
    write_text(ROOT/"logs/PRJNA428535_SUPPORT_COMPLETE",datetime.now(timezone.utc).astimezone().isoformat(timespec="seconds")+"\n")
    print(json.dumps(summary,ensure_ascii=False))


if __name__=="__main__": main()
