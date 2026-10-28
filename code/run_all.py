#!/usr/bin/env python3
"""
run_all.py - regenerate every experiment of "Measurement-Indexed Counting" into a dataset.

Writes, under --out (default ../dataset):
  instances/*.cnf        every formula used, in DIMACS, with generator and seed in comments
  designs/*.json         every measurement design of the toy (groups, assignments, rows, forms)
  results/E1_toy.json    executable model of the OpenAI pipeline (Section 4 of the paper)
  results/E2_components.json   within-group component structure (Section 5)
  results/E3_fibre.json  the fibre loop three ways (Section 6)
  results/E4_constrained.json  semiring component caches vs cardinality encodings (Section 7)
  logs/*.txt             stdout of each experiment
  MANIFEST.json          SHA-256 of every file, environment, timings

All random choices are seeded; every quantity except wall-clock time is deterministic.
--quick skips the two long runs (toy n = 14 with q = 1; encoding at n = 40, b = 20).
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import platform
import random
import sys
import time
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import measurement_toy as T          # noqa: E402
import fibre_loop as FL              # noqa: E402
import constrained_count as CC       # noqa: E402
import numpy as np                   # noqa: E402
import gmpy2                         # noqa: E402


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------
def fr(x):
    """JSON-safe exact number: ints stay ints, Fractions become 'p/q' strings."""
    if isinstance(x, bool):
        return x
    if isinstance(x, int):
        return x
    x = Fraction(x)
    return str(x)


def dimacs(cnf: T.CNF, comments: list[str]) -> str:
    out = [f"c {c}" for c in comments]
    out.append(f"p cnf {cnf.n} {cnf.m}")
    for c in cnf.clauses:
        out.append(" ".join(str(l) for l in c) + " 0")
    return "\n".join(out) + "\n"


def write(path: str, text: str):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as fh:
        fh.write(text)


def sha256(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def hist_digest(hist: dict) -> str:
    items = sorted((list(k), fr(v)) for k, v in hist.items())
    return hashlib.sha256(json.dumps(items).encode()).hexdigest()


class Log:
    def __init__(self, path):
        os.makedirs(os.path.dirname(path), exist_ok=True)
        self.fh = open(path, "w")

    def __call__(self, *a):
        s = " ".join(str(x) for x in a)
        print(s)
        self.fh.write(s + "\n")
        self.fh.flush()

    def close(self):
        self.fh.close()


# ---------------------------------------------------------------------------
# E1: the executable model of the OpenAI pipeline
# ---------------------------------------------------------------------------
def design_json(D: T.ToyDesign, Q: int) -> dict:
    pts = {}
    for t in range(1, Q + 1):
        pt = D.at(t)
        pts[str(t)] = {
            "groups": [[v + 1 for v in S] for S in pt.groups],
            "sigma": [g + 1 for g in pt.sigma],
            "rows": [[{str(v + 1): c for v, c in row.items()} for row in rows_g] for rows_g in pt.rows],
        }
    forms = [[{str(v + 1): [fr(c) for c in cs] for v, cs in D.P[e][k].items()} for k in range(len(D.occ[e]))]
             for e in range(len(D.occ))]
    info = {k: v for k, v in D.info.items() if k not in ("C_rows", "a", "d", "f")}
    if "f" in D.info:
        info["f_coefficients"] = [fr(c) for c in D.info["f"]]
    if "a" in D.info:
        info["variable_schedules"] = [[g + 1 for g in s] for s in D.info["a"]]
        info["clause_schedules"] = [[g + 1 for g in s] for s in D.info["d"]]
    return {"kind": D.kind, "r": D.r, "Q": Q, "degP": D.degP, "M": D.M, "points": pts, "forms": forms, "info": info}


def toy_run(name, cnf, r, kind, deform, q_target, seed, Qaux, Raux, budget, force, log, out_dir):
    rng = random.Random(seed)
    n = cnf.n
    rec = {"name": name, "n": n, "m": cnf.m, "seed": seed, "kind": kind, "r_requested": r,
           "deform": deform, "q_target": q_target, "Qaux": Qaux, "Raux": Raux}
    t0 = time.perf_counter()
    true_count = T.count_solutions(cnf)
    rec["brute_force_count"] = true_count
    rec["t_brute_force_s"] = time.perf_counter() - t0
    if kind == "direct":
        D = T.build_direct_design(cnf, r, deform, rng)
        if D.degP == 0:
            Q, q = 1, 0
        else:
            Q = (D.M * D.degP) // (max(0, q_target) + 1) + 1
            q = (D.M * D.degP) // Q
    elif kind == "collapse":
        D = T.build_collapse_design(cnf, r)
        Q, q = 1, D.M * D.degP
    else:
        D = T.build_aux_design(cnf, r, Qaux, Raux, rng)
        Q, q = D.Q, (D.M * D.degP) // D.Q
    J = q + 1
    degH = D.M * D.degP
    rec.update({"r": D.r, "degP": D.degP, "M": D.M, "degH_bound": degH, "Q": Q, "q": q, "J": J})
    write(os.path.join(out_dir, "designs", f"{name}.json"), json.dumps(design_json(D, Q), indent=1))
    log(f"[{name}] n={n} m={cnf.m} M={D.M} kind={kind} r={D.r} degP={D.degP} Q={Q} q={q} J={J} #F={true_count}")

    counters = T.ToyCounters()
    cache = {}
    points = []
    values = {}
    t26 = time.perf_counter()
    for t in range(1, Q + 1):
        c0 = counters.tuples
        Z = T.measurement_algorithm(D, t, 1, None, frozenset(), counters, cache)
        Ht = T.brute_force_H(cnf, D, t)
        values[t] = Z
        points.append({"t": t, "Z": fr(Z), "H_brute": fr(Ht), "ok": Z == Ht, "tuples": counters.tuples - c0})
        if Z != Ht:
            raise SystemExit(f"identity (2.6) failed in {name} at t={t}")
    rec["identity_2_6"] = {"points": points, "all_ok": all(p["ok"] for p in points),
                           "assignment_evals": counters.assignment_evals, "clause_evals": counters.clause_evals,
                           "t_s": time.perf_counter() - t26}
    pt1 = D.at(1)
    box = 1
    for g in range(len(pt1.groups)):
        Dg = max((sum(abs(c) for c in row.values()) for row in pt1.rows[g]), default=0)
        box *= (2 * Dg + 1) ** len(pt1.rows[g])
    pairs = sum(math.comb(n, s) * (1 << s) for s in range(min(q, n) + 1))
    rec["costs_2_11"] = {"taylor_pairs_per_point": pairs, "tuples_max": max(p["tuples"] for p in points),
                         "box_bound_t1": box, "group_assignments": sum(1 << len(S) for S in pt1.groups), "two_to_n": 1 << n}
    log(f"[{name}] (2.6) OK at {Q} points; tuples {[p['tuples'] for p in points][:6]}...; costs: pairs {pairs} x tuples {rec['costs_2_11']['tuples_max']} x group {rec['costs_2_11']['group_assignments']} vs 2^n {1 << n}")

    H0 = None
    stage = None
    if D.degP == 0:
        H0 = values[1]
        stage = "constant"
    elif q == 0:
        H0 = T.hermite_interpolate({(t, 0): values[t] for t in range(1, Q + 1)}, degH)
        stage = "plain_interpolation"
    elif pairs * Q > budget and not force:
        stage = "skipped_over_budget"
        rec["taylor_stage"] = {"stage": stage, "restricted_sums_needed": pairs * Q, "budget": budget}
    else:
        counters = T.ToyCounters()
        data = {}
        taylor = {}
        t29 = time.perf_counter()
        for t in range(1, Q + 1):
            coeffs = T.taylor_data(D, n, t, q, counters, cache)
            taylor[str(t)] = [fr(c) for c in coeffs]
            for j, hv in enumerate(coeffs):
                data[(t, j)] = hv
        H0 = T.hermite_interpolate(data, degH)
        stage = "taylor_2_9"
        rec["taylor_stage"] = {"stage": stage, "restricted_sums": counters.restricted_sums, "tuples": counters.tuples,
                               "assignment_evals": counters.assignment_evals, "clause_evals": counters.clause_evals,
                               "t_s": time.perf_counter() - t29, "coefficients": taylor, "conditions": len(data)}
    if stage in ("constant", "plain_interpolation"):
        rec["taylor_stage"] = {"stage": stage}
    if H0 is not None:
        rec["H0"] = fr(H0)
        rec["H0_equals_count"] = (H0 == true_count)
        if H0 != true_count:
            raise SystemExit(f"H(0) != #F in {name}")
        log(f"[{name}] H(0) = {H0} = #F ({stage})")
    else:
        log(f"[{name}] Taylor stage skipped ({pairs * Q} restricted sums > budget {budget}); (2.6) certified at every point")
    rec["t_total_s"] = time.perf_counter() - t0
    return rec


def E1(out_dir, log, quick):
    runs = [
        ("toy_direct_n10_q1", 10, 10, 7, dict(r=2, kind="direct", deform=True, q_target=1)),
        ("toy_direct_n12_q1", 12, 12, 3, dict(r=2, kind="direct", deform=True, q_target=1)),
        ("toy_direct_n14_q0", 14, 14, 5, dict(r=2, kind="direct", deform=True, q_target=0)),
        ("toy_direct_n12_nodeform", 12, 12, 3, dict(r=2, kind="direct", deform=False, q_target=0)),
        ("toy_direct_n8_q1", 8, 8, 1, dict(r=2, kind="direct", deform=True, q_target=1)),
        ("toy_aux_n6_Q2R2", 6, 6, 2, dict(r=2, kind="aux", deform=True, q_target=1, Qaux=2, Raux=2)),
        ("toy_aux_n12_Q2R2", 12, 12, 3, dict(r=2, kind="aux", deform=True, q_target=1, Qaux=2, Raux=2)),
        ("toy_aux_n12_Q3R2", 12, 12, 3, dict(r=2, kind="aux", deform=True, q_target=1, Qaux=3, Raux=2)),
        ("toy_aux_n12_r3_Q2R2", 12, 12, 3, dict(r=3, kind="aux", deform=True, q_target=1, Qaux=2, Raux=2)),
        ("toy_aux_n8_Q3R4", 8, 8, 1, dict(r=2, kind="aux", deform=True, q_target=1, Qaux=3, Raux=4)),
        ("toy_collapse_n7", 7, 7, 1, dict(r=2, kind="collapse", deform=True, q_target=1)),
    ]
    if not quick:
        runs.insert(3, ("toy_direct_n14_q1", 14, 14, 5, dict(r=2, kind="direct", deform=True, q_target=1)))
    results = []
    for name, n, m, seed, kw in runs:
        cnf = T.random_bounded_occurrence_3cnf(n, m, random.Random(seed))
        write(os.path.join(out_dir, "instances", f"{name}.cnf"),
              dimacs(cnf, [f"{name}: random bounded-occurrence 3-CNF, measurement_toy.random_bounded_occurrence_3cnf(n={n}, m={m}, seed={seed})",
                           "every variable occurs in at most three clauses; clauses have three distinct variables"]))
        kw = dict(kw)
        kw.setdefault("Qaux", 2)
        kw.setdefault("Raux", 2)
        rec = toy_run(name, cnf, kw["r"], kw["kind"], kw["deform"], kw["q_target"], seed, kw["Qaux"], kw["Raux"],
                      budget=2000, force=False, log=log, out_dir=out_dir)
        results.append(rec)
    # selftest-style check at n = 7 for two seeds: both designs
    checks = []
    for seed in (1, 2):
        cnf = T.random_bounded_occurrence_3cnf(7, 7, random.Random(seed))
        write(os.path.join(out_dir, "instances", f"toy_check_n7_s{seed}.cnf"),
              dimacs(cnf, [f"toy_check_n7_s{seed}: measurement_toy.random_bounded_occurrence_3cnf(n=7, m=7, seed={seed})"]))
        D = T.build_direct_design(cnf, 2, True, random.Random(seed))
        Q = D.M // 2 + 1
        data, cache, ok26 = {}, {}, True
        for t in range(1, Q + 1):
            ok26 &= T.measurement_algorithm(D, t, 1, None, frozenset(), None, cache) == T.brute_force_H(cnf, D, t)
            for j, hv in enumerate(T.taylor_data(D, cnf.n, t, D.M // Q, T.ToyCounters(), cache)):
                data[(t, j)] = hv
        H0 = T.hermite_interpolate(data, D.M)
        A = T.build_aux_design(cnf, 2, 2, 2, random.Random(seed))
        okA = all(T.measurement_algorithm(A, t, 1, None, frozenset(), None, {}) == T.brute_force_H(cnf, A, t) for t in (1, 2))
        checks.append({"seed": seed, "n": 7, "direct_2_6_ok": ok26, "direct_H0": fr(H0), "count": T.count_solutions(cnf),
                       "direct_H0_equals_count": H0 == T.count_solutions(cnf), "aux_2_6_ok": okA})
        log(f"[toy_check_n7_s{seed}] direct (2.6) {ok26}, H0 {H0} == #F {T.count_solutions(cnf)}: {H0 == T.count_solutions(cnf)}; aux (2.6) {okA}")
    return {"runs": results, "checks_n7": checks}


# ---------------------------------------------------------------------------
# E2: within-group component structure under random schedules
# ---------------------------------------------------------------------------
def E2(out_dir, log, quick):
    out = []
    sizes = [1000, 4000] if quick else [1000, 4000, 16000]
    for n in sizes:
        for r in (2, 3, 4, 5, 6, 8, 12):
            for seed in range(3):
                rng = random.Random(seed)
                cnf = T.random_bounded_occurrence_3cnf(n, n, rng)
                occ = [[abs(l) - 1 for l in cl] for cl in cnf.clauses]
                a = [rng.randrange(r) for _ in range(n)]
                d = [rng.randrange(r) for _ in range(len(occ))]
                parent = list(range(n))

                def find(x):
                    while parent[x] != x:
                        parent[x] = parent[parent[x]]
                        x = parent[x]
                    return x
                for e, cl in enumerate(occ):
                    own = [v for v in cl if a[v] == d[e]]
                    for v in own[1:]:
                        parent[find(v)] = find(own[0])
                comp_sizes = {}
                for v in range(n):
                    comp_sizes[find(v)] = comp_sizes.get(find(v), 0) + 1
                sz = sorted(comp_sizes.values(), reverse=True)
                group_sizes = [sum(1 for v in range(n) if a[v] == g) for g in range(r)]
                out.append({"n": n, "m": cnf.m, "r": r, "seed": seed, "largest": sz[0], "second": sz[1] if len(sz) > 1 else 0,
                            "components": len(sz), "singletons": sum(1 for s in sz if s == 1),
                            "mean_group_size": n / r, "largest_over_group": sz[0] / (n / r),
                            "expected_within_group_degree": 3 / r, "group_sizes": group_sizes})
        log(f"[E2] n={n}: largest within-group component by r: " +
            ", ".join(f"r={r}: {min(o['largest'] for o in out if o['n'] == n and o['r'] == r)}..{max(o['largest'] for o in out if o['n'] == n and o['r'] == r)}" for r in (2, 3, 4, 5, 6, 8, 12)))
    return out


# ---------------------------------------------------------------------------
# E3: the fibre loop three ways
# ---------------------------------------------------------------------------
def fibre_run(G: FL.GroupInstance, name, log, out_dir, dict_limit=22, bucket_limit=24, store_hist_k=12):
    rec = {"name": name, "label": G.label, "k": G.k, "m_g": G.m_g, "rows": len(G.rows), "lanes": G.lanes, "radix": G.radix,
           "rows_detail": [{str(v + 1): c for v, c in row.items()} for row in G.rows],
           "prime": FL.PRIME,
           "group": {"variables": [v + 1 for v in G.S],
                     "clauses": [{"signs": list(signs), "forms": [{str(v + 1): c for v, c in f.items()} for f in forms],
                                  "shifts": list(shifts)} for signs, forms, shifts in G.clauses]}}
    comps = G.components()
    sizes = sorted((len(v) for v, _ in comps), reverse=True)
    rec["components"] = {"count": len(comps), "sizes": sizes, "enumerated": sum(1 << s for s in sizes), "two_to_k": 1 << G.k}
    res = {}
    h_bucket = None
    if G.k <= bucket_limit:
        t0 = time.perf_counter(); h_bucket = G.bucket(); res["bucket"] = {"t_s": time.perf_counter() - t0, "attained": len(h_bucket)}
    h_dict = None
    if G.k <= dict_limit:
        t0 = time.perf_counter(); h_dict = G.dict_conv(); res["dict"] = {"t_s": time.perf_counter() - t0, **G.stats_dict}
    t0 = time.perf_counter(); h_packed = G.packed_conv(); res["packed"] = {"t_s": time.perf_counter() - t0, **G.stats_packed}
    rec["methods"] = res
    agree = (h_bucket is None or h_bucket == h_packed) and (h_dict is None or h_dict == h_packed)
    rec["agree"] = agree
    rec["histogram_digest"] = hist_digest(h_packed)
    rec["attained"] = len(h_packed)
    if G.k <= 12 and G.exact_clauses is not None:
        he = G.exact_histogram()
        rec["exact_rational_check"] = (he == h_packed)
        agree = agree and he == h_packed
    if G.k <= store_hist_k:
        rec["histogram"] = sorted((list(k), v) for k, v in h_packed.items())
    if not agree:
        raise SystemExit(f"fibre histograms disagree in {name}")
    tb = res.get("bucket", {}).get("t_s")
    log(f"[{name}] k={G.k} rows={len(G.rows)} lanes={G.lanes} comps={sizes[:6]} enumerated={rec['components']['enumerated']} "
        f"bucket={tb if tb is None else round(tb*1e3,1)}ms dict={round(res['dict']['t_s']*1e3,1) if 'dict' in res else None}ms "
        f"packed={round(res['packed']['t_s']*1e3,1)}ms agree={agree}")
    return rec


def E3(out_dir, log, quick):
    out = {"design_r2": [], "synthetic": []}
    for k in (8, 10, 12, 14, 16, 18, 20):
        n = 2 * k
        G = FL.GroupInstance.from_design(n, 2, 1)
        cnf = T.random_bounded_occurrence_3cnf(n, n, random.Random(1))
        write(os.path.join(out_dir, "instances", f"fibre_design_r2_k{k}.cnf"),
              dimacs(cnf, [f"fibre_design_r2_k{k}: measurement_toy.random_bounded_occurrence_3cnf(n={n}, m={n}, seed=1); aux design r=2 Q=2 R=2 built from the same seed stream"]))
        out["design_r2"].append(fibre_run(G, f"fibre_design_r2_k{k}", log, out_dir))
    synth = [(4, 8, (16, 18, 20, 22, 24)), (3, 8, (16, 18, 20, 22)), (4, 10, (20, 22, 24))]
    for r, rows, ks in synth:
        for k in ks:
            if quick and k >= 24:
                continue
            n = r * k
            G = FL.GroupInstance.synthetic(n, r, rows, 1, 3)
            cnf = T.random_bounded_occurrence_3cnf(n, n, random.Random(1))
            write(os.path.join(out_dir, "instances", f"fibre_synth_r{r}_rows{rows}_k{k}.cnf"),
                  dimacs(cnf, [f"fibre_synth_r{r}_rows{rows}_k{k}: measurement_toy.random_bounded_occurrence_3cnf(n={n}, m={n}, seed=1); group 0 = {k} random variables, random clause assignment, {rows} synthetic +-1 rows of support 3"]))
            rec = fibre_run(G, f"fibre_synth_r{r}_rows{rows}_k{k}", log, out_dir, dict_limit=24 if rows == 10 else 22)
            rec.update({"r": r, "synthetic_rows": rows})
            out["synthetic"].append(rec)
    return out


# ---------------------------------------------------------------------------
# E4: semiring component caches vs cardinality encodings
# ---------------------------------------------------------------------------
def cc_stats(C: CC.Counter) -> dict:
    s = C.stats
    return {"nodes": s.nodes, "cache_entries": s.cache_misses, "cache_hits": s.cache_hits, "conflicts": s.conflicts,
            "max_component": s.max_component, "hist_entries_cached": s.hist_entries_cached, "convolutions": s.convolutions}


def E4(out_dir, log, quick):
    out = {"grids": [], "multirow": [], "random": []}
    for W, H in ((4, 4), (5, 5), (6, 6), (7, 7), (8, 8)):
        n, clauses = CC.grid_independent_sets(W, H)
        write(os.path.join(out_dir, "instances", f"grid_{W}x{H}.cnf"),
              dimacs(T.CNF(n, clauses), [f"grid_{W}x{H}: independent sets of the {W}x{H} grid graph; vertices row-major, one binary clause per edge"]))
        ones = [[1] * n]
        C = CC.Counter(n, clauses, rows=ones, order="static")
        t0 = time.perf_counter(); hist = C.run(); dt = time.perf_counter() - t0
        total = sum(hist.values())
        rec = {"grid": [W, H], "n": n, "m": len(clauses), "histogram": sorted((k[0], v) for k, v in hist.items()),
               "total": total, "oeis_A006506": CC.OEIS_A006506.get(W) if W == H else None,
               "total_equals_oeis": (total == CC.OEIS_A006506.get(W)) if W == H else None,
               "t_s": dt, **cc_stats(C), "encodings": []}
        if n <= 20:
            rec["brute_force_ok"] = (CC.brute_force(n, clauses, C.w0, C.w1, ones) == hist)
        mode = max(hist, key=lambda y: hist[y])[0]
        rec["mode"] = mode
        log(f"[grid {W}x{H}] n={n} sizes={len(hist)} total={total} oeis={rec['total_equals_oeis']} nodes={C.stats.nodes} t={dt*1e3:.1f}ms")
        if (W, H) in ((4, 4), (5, 5), (6, 6)):
            for b in sorted({mode, mode // 2, 1}):
                n2, cl2, aux = CC.encode_cardinality(n, range(1, n + 1), b, clauses)
                C2 = CC.Counter(n2, cl2, order="static")
                t0 = time.perf_counter(); h2 = C2.run(); dt2 = time.perf_counter() - t0
                val = h2.get((), 0)
                rec["encodings"].append({"b": b, "aux": aux, "clauses": len(cl2), "count": val, "equals_histogram": val == hist.get((b,), 0),
                                         "t_s": dt2, "ratio_to_histogram": dt2 / dt, **cc_stats(C2)})
                log(f"   encoding b={b}: aux={aux} nodes={C2.stats.nodes} maxcomp={C2.stats.max_component} t={dt2*1e3:.1f}ms ratio={dt2/dt:.2f} ok={val == hist.get((b,), 0)}")
        out["grids"].append(rec)
    # multirow at 6x6
    rng = random.Random(1)
    n, clauses = CC.grid_independent_sets(6, 6)
    pb = [[rng.randint(1, 3) for _ in range(n)]]
    C = CC.Counter(n, clauses, rows=pb, order="static")
    t0 = time.perf_counter(); hist = C.run(); dt = time.perf_counter() - t0
    out["multirow"].append({"grid": [6, 6], "kind": "pseudo_boolean", "row": pb[0], "histogram": sorted((k[0], v) for k, v in hist.items()),
                            "t_s": dt, **cc_stats(C)})
    left = [1 if (v % 6) < 3 else 0 for v in range(n)]
    right = [1 - l for l in left]
    C = CC.Counter(n, clauses, rows=[left, right], order="static")
    t0 = time.perf_counter(); hist2 = C.run(); dt2 = time.perf_counter() - t0
    out["multirow"].append({"grid": [6, 6], "kind": "two_rows", "rows": [left, right], "histogram": sorted((list(k), v) for k, v in hist2.items()),
                            "t_s": dt2, **cc_stats(C)})
    log(f"[6x6 multirow] pseudo-Boolean {len(hist)} values {dt*1e3:.1f}ms; two rows {len(hist2)} pairs {dt2*1e3:.1f}ms")
    # random 3-CNF
    for n, m, weighted, bs in ((20, 50, False, [5, 10]), (30, 75, True, [7, 15]), (40, 100, True, [10, 20])):
        rng = random.Random(1)
        clauses = CC.random_3cnf(n, m, rng)
        w0 = [1] * n
        w1 = [rng.randint(1, 3) for _ in range(n)] if weighted else [1] * n
        write(os.path.join(out_dir, "instances", f"random3cnf_n{n}_m{m}_s1.cnf"),
              dimacs(T.CNF(n, clauses), [f"random3cnf_n{n}_m{m}_s1: constrained_count.random_3cnf(n={n}, m={m}, seed=1)",
                                         f"literal weights w1 = {w1} (w0 = 1)" if weighted else "unweighted"]))
        ones = [[1] * n]
        C = CC.Counter(n, clauses, w0, w1, rows=ones, order="dynamic")
        t0 = time.perf_counter(); hist = C.run(); dt = time.perf_counter() - t0
        C0 = CC.Counter(n, clauses, w0, w1, order="dynamic")
        t0 = time.perf_counter(); plain = C0.run(); dt0 = time.perf_counter() - t0
        rec = {"n": n, "m": m, "weighted": weighted, "w1": w1, "histogram": sorted((k[0], v) for k, v in hist.items()),
               "total": sum(hist.values()), "plain_count": plain.get((), 0), "plain_equals_total": plain.get((), 0) == sum(hist.values()),
               "t_hist_s": dt, "t_plain_s": dt0, "hist_stats": cc_stats(C), "plain_stats": cc_stats(C0), "encodings": []}
        if n <= 20:
            rec["brute_force_ok"] = (CC.brute_force(n, clauses, w0, w1, ones) == hist)
        log(f"[random n={n}] hist {len(hist)} entries t={dt*1e3:.1f}ms nodes={C.stats.nodes}; plain t={dt0*1e3:.1f}ms")
        for b in bs:
            if quick and n == 40 and b == 20:
                continue
            n2, cl2, aux = CC.encode_cardinality(n, range(1, n + 1), b, clauses)
            C2 = CC.Counter(n2, cl2, w0 + [1] * (n2 - n), w1 + [1] * (n2 - n), order="dynamic")
            t0 = time.perf_counter(); h2 = C2.run(); dt2 = time.perf_counter() - t0
            val = h2.get((), 0)
            rec["encodings"].append({"b": b, "aux": aux, "clauses": len(cl2), "count": val, "equals_histogram": val == hist.get((b,), 0),
                                     "t_s": dt2, "ratio_to_histogram": dt2 / dt, **cc_stats(C2)})
            log(f"   encoding b={b}: aux={aux} nodes={C2.stats.nodes} maxcomp={C2.stats.max_component} t={dt2*1e3:.1f}ms ratio={dt2/dt:.1f} ok={val == hist.get((b,), 0)}")
        out["random"].append(rec)
    return out


# ---------------------------------------------------------------------------
def environment() -> dict:
    return {"python": sys.version, "platform": platform.platform(), "machine": platform.machine(),
            "cpu_count": os.cpu_count(), "numpy": np.__version__, "gmpy2": gmpy2.version(),
            "started_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=os.path.join(HERE, "..", "dataset"))
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--only", default="", help="comma-separated subset of E1,E2,E3,E4")
    ap.add_argument("--manifest-only", action="store_true", help="only (re)write MANIFEST.json for an existing dataset")
    args = ap.parse_args()
    out_dir = os.path.abspath(args.out)
    os.makedirs(out_dir, exist_ok=True)
    if args.manifest_only:
        old = json.load(open(os.path.join(out_dir, "MANIFEST.json"))) if os.path.exists(os.path.join(out_dir, "MANIFEST.json")) else {}
        write_manifest(out_dir, old.get("environment", environment()), old.get("quick", False), old.get("timings_s", {}))
        print("manifest rewritten")
        return
    env = environment()
    timings = {}
    only = set(args.only.split(",")) if args.only else {"E1", "E2", "E3", "E4"}
    for tag, fn in (("E1", E1), ("E2", E2), ("E3", E3), ("E4", E4)):
        if tag not in only:
            continue
        log = Log(os.path.join(out_dir, "logs", f"{tag}.txt"))
        t0 = time.perf_counter()
        res = fn(out_dir, log, args.quick)
        timings[tag] = time.perf_counter() - t0
        log(f"[{tag}] done in {timings[tag]:.1f} s")
        log.close()
        names = {"E1": "E1_toy", "E2": "E2_components", "E3": "E3_fibre", "E4": "E4_constrained"}
        write(os.path.join(out_dir, "results", f"{names[tag]}.json"),
              json.dumps({"experiment": tag, "environment": env, "quick": args.quick, "data": res}, indent=1))
    write_manifest(out_dir, env, args.quick, timings)


def write_manifest(out_dir, env, quick, timings):
    # copy the dataset's README and LICENSE in before hashing
    import shutil
    meta = os.path.join(HERE, "..", "dataset_meta")
    if os.path.isdir(meta):
        for f in os.listdir(meta):
            shutil.copy(os.path.join(meta, f), os.path.join(out_dir, f))
    manifest = {"environment": env, "quick": quick, "timings_s": timings, "files": {}}
    for root, _, files in os.walk(out_dir):
        for f in sorted(files):
            if f == "MANIFEST.json":
                continue
            p = os.path.join(root, f)
            manifest["files"][os.path.relpath(p, out_dir)] = {"sha256": sha256(p), "bytes": os.path.getsize(p)}
    write(os.path.join(out_dir, "MANIFEST.json"), json.dumps(manifest, indent=1))
    print(f"dataset written to {out_dir}: {len(manifest['files'])} files; timings {timings}")


if __name__ == "__main__":
    main()
