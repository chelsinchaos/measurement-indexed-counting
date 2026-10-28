#!/usr/bin/env python3
"""
verify_all.py - one command for an auditor.

1. Regenerates the experiments with run_all.py --quick into a fresh directory and compares every
   deterministic quantity (counts, histograms, digests, node counts, design points, Taylor data)
   with the released dataset; wall-clock times are ignored.
2. Runs independent_check.py (different algorithms) against the released dataset.
3. Recomputes the SHA-256 manifest of the released dataset and reports mismatches.

Exit code 0 only if every comparison and check passes.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
TIMING_KEYS = {"t_s", "t_total_s", "t_brute_force_s", "t_hist_s", "t_plain_s", "ratio_to_histogram",
               "t_kernel", "t_mul", "t_reduce", "environment", "quick", "timings_s", "started_utc", "finished_utc", "label"}


def strip(o):
    if isinstance(o, dict):
        return {k: strip(v) for k, v in o.items() if k not in TIMING_KEYS}
    if isinstance(o, list):
        return [strip(v) for v in o]
    return o


def load(ds, name):
    return json.load(open(os.path.join(ds, "results", name)))["data"]


def compare(a, b, path, problems):
    a, b = strip(a), strip(b)
    if a != b:
        if isinstance(a, dict) and isinstance(b, dict):
            for k in sorted(set(a) | set(b)):
                if k not in a or k not in b:
                    problems.append(f"{path}.{k}: present in only one side")
                elif a[k] != b[k]:
                    compare(a[k], b[k], f"{path}.{k}", problems)
        elif isinstance(a, list) and isinstance(b, list) and len(a) == len(b):
            for i, (x, y) in enumerate(zip(a, b)):
                if x != y:
                    compare(x, y, f"{path}[{i}]", problems)
        else:
            problems.append(f"{path}: {str(a)[:80]} != {str(b)[:80]}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", default=os.path.join(HERE, "..", "dataset"))
    ap.add_argument("--skip-regeneration", action="store_true")
    ap.add_argument("--skip-independent", action="store_true")
    args = ap.parse_args()
    ds = os.path.abspath(args.dataset)
    ok = True

    # 3. manifest
    man = json.load(open(os.path.join(ds, "MANIFEST.json")))
    bad = 0
    for rel, info in man["files"].items():
        p = os.path.join(ds, rel)
        if not os.path.exists(p):
            print(f"MANIFEST: missing {rel}"); bad += 1; continue
        h = hashlib.sha256(open(p, "rb").read()).hexdigest()
        if h != info["sha256"]:
            print(f"MANIFEST: sha256 mismatch {rel}"); bad += 1
    print(f"manifest: {len(man['files'])} files, {bad} problems")
    ok &= bad == 0

    # 1. regeneration
    if not args.skip_regeneration:
        tmp = tempfile.mkdtemp(prefix="mic_regen_")
        print(f"regenerating with --quick into {tmp} ...")
        r = subprocess.run([sys.executable, os.path.join(HERE, "run_all.py"), "--quick", "--out", tmp],
                           capture_output=True, text=True)
        if r.returncode != 0:
            print("regeneration FAILED:\n" + r.stdout[-3000:] + r.stderr[-3000:]); ok = False
        else:
            problems = []
            # E1 by run name (the quick run lacks toy_direct_n14_q1)
            A = {x["name"]: x for x in load(ds, "E1_toy.json")["runs"]}
            B = {x["name"]: x for x in load(tmp, "E1_toy.json")["runs"]}
            for name in sorted(set(A) & set(B)):
                compare(A[name], B[name], f"E1.{name}", problems)
            compare(load(ds, "E1_toy.json")["checks_n7"], load(tmp, "E1_toy.json")["checks_n7"], "E1.checks_n7", problems)
            # E2 by (n, r, seed) for n <= 4000
            A = {(x["n"], x["r"], x["seed"]): x for x in load(ds, "E2_components.json") if x["n"] <= 4000}
            B = {(x["n"], x["r"], x["seed"]): x for x in load(tmp, "E2_components.json")}
            for key in sorted(set(A) & set(B)):
                compare(A[key], B[key], f"E2.{key}", problems)
            # E3 by name
            A = {x["name"]: x for x in load(ds, "E3_fibre.json")["design_r2"] + load(ds, "E3_fibre.json")["synthetic"]}
            B = {x["name"]: x for x in load(tmp, "E3_fibre.json")["design_r2"] + load(tmp, "E3_fibre.json")["synthetic"]}
            for name in sorted(set(A) & set(B)):
                compare(A[name], B[name], f"E3.{name}", problems)
            # E4
            A, B = load(ds, "E4_constrained.json"), load(tmp, "E4_constrained.json")
            for ga, gb in zip(A["grids"], B["grids"]):
                compare(ga, gb, f"E4.grid{ga['grid']}", problems)
            for ma, mb in zip(A["multirow"], B["multirow"]):
                compare(ma, mb, f"E4.multirow.{ma['kind']}", problems)
            for ra, rb in zip(A["random"], B["random"]):
                ea = {e["b"]: e for e in ra["encodings"]}
                eb = {e["b"]: e for e in rb["encodings"]}
                ra2 = {k: v for k, v in ra.items() if k != "encodings"}
                rb2 = {k: v for k, v in rb.items() if k != "encodings"}
                compare(ra2, rb2, f"E4.random{ra['n']}", problems)
                for b in sorted(set(ea) & set(eb)):
                    compare(ea[b], eb[b], f"E4.random{ra['n']}.encoding{b}", problems)
            print(f"regeneration: {len(problems)} differences in deterministic quantities")
            for p in problems[:50]:
                print("  " + p)
            ok &= not problems

    # 2. independent check
    if not args.skip_independent:
        r = subprocess.run([sys.executable, os.path.join(HERE, "independent_check.py"), "--dataset", ds],
                           capture_output=True, text=True, cwd=tempfile.mkdtemp(prefix="mic_indep_"))
        print(r.stdout[-6000:])
        if r.returncode != 0:
            print("independent check FAILED"); ok = False
    print("VERIFY_ALL: " + ("PASS" if ok else "FAIL"))
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
