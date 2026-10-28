#!/usr/bin/env python3
"""
independent_check.py - verification run 2: recompute the results of the dataset by DIFFERENT
algorithms, written separately from the code that produced them, and compare.

  E1 (toy)        H(T) = sum_x prod_e Phi_e(<P_{e,v}(T), x>) is computed as an explicit polynomial
                  in T by symbolic brute force over all 2^n assignments (no measurements, no
                  designs, no interpolation). Its values at t = 1..Q must equal the recorded
                  measurement-algorithm values Z; its Taylor coefficients at t must equal the
                  recorded ones; its constant term must equal an independent model count
                  (exhaustive DPLL). Design condition (4) of Definition 2.1 is re-tested with a
                  fresh rank computation; conditions (1) and (3) are re-measured.
  E2 (components) recomputed with a breadth-first search instead of union-find.
  E3 (fibre loop) every histogram with k <= 16 recomputed mod p by plain enumeration in pure
                  Python from the stored clauses, forms, shifts and rows; compared with the
                  stored histogram (k <= 12) or digest.
  E4 (constrained) grid histograms recomputed by a transfer-matrix (sliding-window) dynamic
                  programme that never decomposes or caches components; random 3-CNF histograms
                  recomputed by exhaustive DPLL enumeration that visits every model; grid totals
                  compared with OEIS A006506; encoded counts compared with histogram entries.

Prints one line per check and a PASS/FAIL summary; writes results/independent_check.json.
Standard library only.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
import time
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))


def frac(s):
    return Fraction(s) if isinstance(s, str) else Fraction(int(s))


# ---------------------------------------------------------------------------
# polynomials in T with Fraction coefficients (independent, minimal)
# ---------------------------------------------------------------------------
def padd(a, b):
    n = max(len(a), len(b))
    return [(a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0) for i in range(n)]


def pmul(a, b):
    out = [Fraction(0)] * (len(a) + len(b) - 1) if a and b else []
    for i, x in enumerate(a):
        if x == 0:
            continue
        for j, y in enumerate(b):
            out[i + j] += x * y
    return out


def pscale(a, c):
    return [x * c for x in a]


def peval(a, t):
    return sum((c * Fraction(t) ** i for i, c in enumerate(a)), Fraction(0))


def pshift(a, t):
    """Coefficients of p(t + z) in z."""
    from math import comb
    out = [Fraction(0)] * len(a)
    for j, c in enumerate(a):
        for i in range(j + 1):
            out[i] += c * comb(j, i) * Fraction(t) ** (j - i)
    return out


def trim(a):
    a = list(a)
    while a and a[-1] == 0:
        a.pop()
    return a


# ---------------------------------------------------------------------------
# DIMACS and an exhaustive DPLL enumerator (visits every model; no cache, no decomposition)
# ---------------------------------------------------------------------------
def read_dimacs(path):
    n = 0
    clauses = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line[0] == "c":
                continue
            if line.startswith("p"):
                n = int(line.split()[2])
                continue
            lits = [int(t) for t in line.split()]
            if lits and lits[-1] == 0:
                lits = lits[:-1]
            if lits:
                clauses.append(tuple(lits))
    return n, clauses


def dpll_enumerate(n, clauses, weight_fn, key_fn):
    """Exhaustive enumeration of all models by DPLL with unit propagation. For every model x
    (a dict var -> 0/1) adds weight_fn(x) to hist[key_fn(x)]. Returns hist."""
    sys.setrecursionlimit(10000)
    hist = {}
    assign = {}

    def propagate(cls):
        changed = True
        cur = cls
        while changed:
            changed = False
            nxt = []
            for cl in cur:
                sat = False
                rest = []
                for l in cl:
                    v = abs(l)
                    if v in assign:
                        if assign[v] == (l > 0):
                            sat = True
                            break
                    else:
                        rest.append(l)
                if sat:
                    continue
                if not rest:
                    return None
                if len(rest) == 1:
                    assign[abs(rest[0])] = rest[0] > 0
                    changed = True
                    continue
                nxt.append(tuple(rest))
            cur = nxt
        return cur

    def rec(cls, free):
        saved = dict(assign)
        cur = propagate(cls)
        if cur is None:
            assign.clear(); assign.update(saved)
            return
        unassigned = [v for v in free if v not in assign]
        if not cur:
            # every remaining variable is free: enumerate them all (they are independent)
            def fill(i):
                if i == len(unassigned):
                    x = {v: (1 if assign[v] else 0) for v in range(1, n + 1)}
                    k = key_fn(x)
                    hist[k] = hist.get(k, 0) + weight_fn(x)
                    return
                for val in (True, False):
                    assign[unassigned[i]] = val
                    fill(i + 1)
                    del assign[unassigned[i]]
            fill(0)
            assign.clear(); assign.update(saved)
            return
        v = abs(cur[0][0])
        for val in (True, False):
            assign[v] = val
            rec(cur, unassigned)
            del assign[v]
        assign.clear(); assign.update(saved)

    rec(list(clauses), list(range(1, n + 1)))
    return hist


# ---------------------------------------------------------------------------
# independent rank test over Q
# ---------------------------------------------------------------------------
def rank_Q(rows):
    """Rank of a list of dict vectors over Q by elimination (fresh implementation)."""
    cols = sorted({k for r in rows for k in r})
    M = [[Fraction(r.get(c, 0)) for c in cols] for r in rows]
    rank = 0
    for c in range(len(cols)):
        piv = next((i for i in range(rank, len(M)) if M[i][c] != 0), None)
        if piv is None:
            continue
        M[rank], M[piv] = M[piv], M[rank]
        pv = M[rank][c]
        M[rank] = [x / pv for x in M[rank]]
        for i in range(len(M)):
            if i != rank and M[i][c] != 0:
                f = M[i][c]
                M[i] = [a - f * b for a, b in zip(M[i], M[rank])]
        rank += 1
    return rank


# ---------------------------------------------------------------------------
# E1
# ---------------------------------------------------------------------------
def check_E1(ds, report):
    data = json.load(open(os.path.join(ds, "results", "E1_toy.json")))["data"]
    for run in data["runs"]:
        name = run["name"]
        n, clauses = read_dimacs(os.path.join(ds, "instances", f"{name}.cnf"))
        design = json.load(open(os.path.join(ds, "designs", f"{name}.json")))
        forms = design["forms"]                       # forms[e][k]: {var: [coeffs]}
        occ = [[(abs(l), l > 0) for l in cl] for cl in clauses]
        # symbolic brute force: H(T) as a polynomial
        t0 = time.time()
        H = [Fraction(0)]
        for bits in range(1 << n):
            x = [(bits >> i) & 1 for i in range(n)]
            prod = [Fraction(1)]
            for e, cl in enumerate(occ):
                # Phi_e = 1 - prod_k (1 - L_k), L_k = input (positive) or 1 - input
                pk = [Fraction(1)]
                for k, (v, pos) in enumerate(cl):
                    inp = [Fraction(0)]
                    for var, cs in forms[e][k].items():
                        if x[int(var) - 1]:
                            inp = padd(inp, [frac(c) for c in cs])
                    L = inp if pos else padd([Fraction(1)], pscale(inp, -1))
                    pk = pmul(pk, padd([Fraction(1)], pscale(L, -1)))
                phi = padd([Fraction(1)], pscale(pk, -1))
                prod = pmul(prod, phi)
                if not trim(prod):
                    break
            H = padd(H, prod)
        H = trim(H)
        dt = time.time() - t0
        degH = len(H) - 1
        # (a) values at the design points
        ok_vals = all(peval(H, p["t"]) == frac(p["Z"]) for p in run["identity_2_6"]["points"])
        report(f"E1 {name}: symbolic H(T) of degree {degH} (bound {run['degH_bound']}) in {dt:.1f}s; "
               f"H(t) equals recorded Z at all {len(run['identity_2_6']['points'])} points: {ok_vals}", ok_vals and degH <= run["degH_bound"])
        # (b) Taylor coefficients
        ts = run.get("taylor_stage", {})
        if ts.get("stage") == "taylor_2_9":
            ok_t = True
            for t, coeffs in ts["coefficients"].items():
                sh = pshift(H, int(t))
                sh = sh + [Fraction(0)] * (len(coeffs) - len(sh))
                ok_t &= all(sh[j] == frac(c) for j, c in enumerate(coeffs))
            report(f"E1 {name}: recorded Taylor coefficients at {len(ts['coefficients'])} points equal the symbolic expansion: {ok_t}", ok_t)
        # (c) constant term = independent model count
        cnt = dpll_enumerate(n, clauses, lambda x: 1, lambda x: ())
        cnt = cnt.get((), 0)
        ok_c = (H[0] if H else 0) == cnt == run["brute_force_count"] and (("H0" not in run) or frac(run["H0"]) == cnt)
        report(f"E1 {name}: H(0) = {H[0] if H else 0}; exhaustive DPLL count = {cnt}; recorded #F = {run['brute_force_count']}"
               + (f"; recorded H(0) = {run['H0']}" if "H0" in run else "") + f": {ok_c}", ok_c)
        # (d) design conditions at every point with a fresh rank test
        ok_d = True
        max_abs_sum = 0
        for t, pt in design["points"].items():
            t = int(t)
            groups = pt["groups"]
            group_of = {}
            for g, S in enumerate(groups):
                for v in S:
                    group_of[v] = g
            # (1) group sizes
            if any(len(S) > 2 * n / len(groups) for S in groups):
                ok_d = False
            for g, rows in enumerate(pt["rows"]):
                for row in rows:
                    max_abs_sum = max(max_abs_sum, sum(abs(c) for c in row.values()))
            for e, cl in enumerate(occ):
                g = pt["sigma"][e] - 1
                for k, (v, pos) in enumerate(cl):
                    vec = {int(var): peval([frac(c) for c in cs], t) for var, cs in forms[e][k].items()}
                    # P(0) = e_v
                    at0 = {int(var): frac(cs[0]) for var, cs in forms[e][k].items() if frac(cs[0]) != 0}
                    if at0 != {v: Fraction(1)}:
                        ok_d = False
                    for h in range(len(groups)):
                        if h == g:
                            continue
                        rest = {var: c for var, c in vec.items() if group_of[var] == h and c != 0}
                        if not rest:
                            continue
                        rows = [{int(var): c for var, c in row.items()} for row in pt["rows"][h]]
                        if rank_Q(rows + [rest]) != rank_Q(rows):
                            ok_d = False
        report(f"E1 {name}: design conditions (1), (3) [max |row| sum {max_abs_sum}], (4) and P(0) = e_v re-verified at every point: {ok_d}", ok_d)


# ---------------------------------------------------------------------------
# E2
# ---------------------------------------------------------------------------
def check_E2(ds, report):
    import random
    sys.path.insert(0, HERE)
    import measurement_toy as T
    data = json.load(open(os.path.join(ds, "results", "E2_components.json")))["data"]
    ok_all = True
    checked = 0
    for rec in data:
        if rec["n"] > 4000:
            continue
        rng = random.Random(rec["seed"])
        cnf = T.random_bounded_occurrence_3cnf(rec["n"], rec["n"], rng)
        a = [rng.randrange(rec["r"]) for _ in range(rec["n"])]
        d = [rng.randrange(rec["r"]) for _ in range(cnf.m)]
        adj = {v: set() for v in range(rec["n"])}
        for e, cl in enumerate(cnf.clauses):
            own = [abs(l) - 1 for l in cl if a[abs(l) - 1] == d[e]]
            for u in own:
                for w in own:
                    if u != w:
                        adj[u].add(w)
        seen = set()
        largest = 0
        comps = 0
        for s in range(rec["n"]):
            if s in seen:
                continue
            comps += 1
            stack = [s]
            seen.add(s)
            size = 0
            while stack:
                u = stack.pop()
                size += 1
                for w in adj[u]:
                    if w not in seen:
                        seen.add(w)
                        stack.append(w)
            largest = max(largest, size)
        ok = largest == rec["largest"] and comps == rec["components"]
        ok_all &= ok
        checked += 1
    report(f"E2: BFS recomputation of largest component and component count agrees on {checked} (n, r, seed) records: {ok_all}", ok_all)


# ---------------------------------------------------------------------------
# E3
# ---------------------------------------------------------------------------
def check_E3(ds, report):
    data = json.load(open(os.path.join(ds, "results", "E3_fibre.json")))["data"]
    recs = data["design_r2"] + data["synthetic"]
    for rec in recs:
        k = rec["k"]
        if k > 16:
            continue
        p = rec["prime"]
        S = rec["group"]["variables"]
        pos = {v: i for i, v in enumerate(S)}
        rows = [{int(v): c for v, c in row.items()} for row in rec["rows_detail"]]
        clauses = [(c["signs"], [{int(v): cc for v, cc in f.items()} for f in c["forms"]], c["shifts"]) for c in rec["group"]["clauses"]]
        hist = {}
        t0 = time.time()
        for bits in range(1 << k):
            w = 1
            for signs, forms, shifts in clauses:
                prod = 1
                for s, f, sh in zip(signs, forms, shifts):
                    X = (sh + sum(c * ((bits >> pos[v]) & 1) for v, c in f.items())) % p
                    L = X if s else (1 - X) % p
                    prod = prod * ((1 - L) % p) % p
                w = w * ((1 - prod) % p) % p
                if w == 0:
                    break
            if w == 0:
                continue
            y = tuple(sum(c * ((bits >> pos[v]) & 1) for v, c in row.items()) for row in rows)
            hist[y] = (hist.get(y, 0) + w) % p
        hist = {y: v for y, v in hist.items() if v}
        items = sorted((list(y), v) for y, v in hist.items())
        digest = hashlib.sha256(json.dumps(items).encode()).hexdigest()
        ok = digest == rec["histogram_digest"] and len(hist) == rec["attained"]
        if "histogram" in rec:
            ok &= items == [(list(y), v) for y, v in rec["histogram"]]
        report(f"E3 {rec['name']}: plain enumeration of 2^{k} assignments mod p reproduces the stored histogram ({len(hist)} entries, {time.time()-t0:.1f}s): {ok}", ok)


# ---------------------------------------------------------------------------
# E4
# ---------------------------------------------------------------------------
OEIS_A006506 = {1: 2, 2: 7, 3: 63, 4: 1234, 5: 55447, 6: 5598861, 7: 1280128950, 8: 660647962955}


def grid_transfer_matrix(W, H, rows_meas):
    """Histogram of independent sets of the W x H grid over the measurement rows (list of
    coefficient lists indexed by vertex), by a sliding-window DP over cells in row-major order.
    State: tuple of the last W cell values. No components, no caching of subformulas."""
    c = len(rows_meas)
    zero = tuple(0 for _ in range(c))
    states = {tuple([0] * W): {zero: 1}}
    for r in range(H):
        for col in range(W):
            v = r * W + col
            vec = tuple(int(row[v]) for row in rows_meas)
            new = {}
            for win, hist in states.items():
                up = win[0]                       # the cell W positions back (above), window oldest-first
                left = win[-1] if col > 0 else 0  # the previous cell, if in the same row
                for val in (0, 1):
                    if val == 1 and (up == 1 or left == 1):
                        continue
                    nwin = win[1:] + (val,)
                    tgt = new.setdefault(nwin, {})
                    for y, w in hist.items():
                        ny = tuple(a + b * val for a, b in zip(y, vec)) if val else y
                        tgt[ny] = tgt.get(ny, 0) + w
            states = new
    total = {}
    for hist in states.values():
        for y, w in hist.items():
            total[y] = total.get(y, 0) + w
    return total


def check_E4(ds, report):
    data = json.load(open(os.path.join(ds, "results", "E4_constrained.json")))["data"]
    for rec in data["grids"]:
        W, H = rec["grid"]
        n = W * H
        t0 = time.time()
        hist = grid_transfer_matrix(W, H, [[1] * n])
        got = sorted((y[0], v) for y, v in hist.items())
        ok = got == [tuple(x) for x in rec["histogram"]] or got == [(a, b) for a, b in rec["histogram"]]
        total = sum(hist.values())
        ok_oeis = (total == OEIS_A006506.get(W)) if W == H else True
        report(f"E4 grid {W}x{H}: transfer-matrix DP reproduces the size histogram ({len(got)} sizes, {time.time()-t0:.1f}s): {ok}; total {total} = OEIS: {ok_oeis}", ok and ok_oeis)
        for enc in rec["encodings"]:
            okb = enc["count"] == dict(got).get(enc["b"], 0)
            report(f"E4 grid {W}x{H}: encoded count at b = {enc['b']} ({enc['count']}) equals transfer-matrix entry: {okb}", okb)
    for rec in data["multirow"]:
        W, H = rec["grid"]
        n = W * H
        rows = [rec["row"]] if rec["kind"] == "pseudo_boolean" else rec["rows"]
        hist = grid_transfer_matrix(W, H, rows)
        got = sorted((list(y), v) for y, v in hist.items())
        stored = sorted(([y] if isinstance(y, int) else list(y), v) for y, v in rec["histogram"])
        ok = got == stored
        report(f"E4 6x6 {rec['kind']}: transfer-matrix DP reproduces the {len(got)}-entry histogram: {ok}", ok)
    for rec in data["random"]:
        n, m = rec["n"], rec["m"]
        _, clauses = read_dimacs(os.path.join(ds, "instances", f"random3cnf_n{n}_m{m}_s1.cnf"))
        w1 = rec["w1"]
        t0 = time.time()

        def weight(x, w1=w1):
            w = 1
            for v in range(1, n + 1):
                if x[v]:
                    w *= w1[v - 1]
            return w
        hist = dpll_enumerate(n, clauses, weight, lambda x: (sum(x.values()),))
        got = sorted((y[0], v) for y, v in hist.items())
        ok = got == [(a, b) for a, b in rec["histogram"]]
        report(f"E4 random n={n} m={m}: exhaustive DPLL enumeration of every model reproduces the cardinality histogram "
               f"({len(got)} entries, total {sum(hist.values())}, {time.time()-t0:.1f}s): {ok}", ok)
        for enc in rec["encodings"]:
            okb = enc["count"] == dict(got).get(enc["b"], 0)
            report(f"E4 random n={n}: encoded count at b = {enc['b']} ({enc['count']}) equals the enumeration: {okb}", okb)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", default=os.path.join(HERE, "..", "dataset"))
    ap.add_argument("--only", default="")
    ap.add_argument("--write", action="store_true", help="write results/independent_check.json into the dataset (release generation only)")
    args = ap.parse_args()
    ds = os.path.abspath(args.dataset)
    results = []

    def report(msg, ok):
        results.append({"check": msg, "ok": bool(ok)})
        print(("PASS " if ok else "FAIL ") + msg)
    only = set(args.only.split(",")) if args.only else {"E1", "E2", "E3", "E4"}
    for tag, fn in (("E1", check_E1), ("E2", check_E2), ("E3", check_E3), ("E4", check_E4)):
        if tag in only:
            fn(ds, report)
    n_ok = sum(1 for r in results if r["ok"])
    summary = {"checks": len(results), "passed": n_ok, "failed": len(results) - n_ok, "results": results,
               "finished_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}
    out = os.path.join(ds, "results", "independent_check.json") if args.write else os.path.abspath("independent_check_rerun.json")
    with open(out, "w") as fh:
        json.dump(summary, fh, indent=1)
    print(f"wrote {out}")
    print(f"\n{n_ok} of {len(results)} independent checks passed")
    sys.exit(0 if n_ok == len(results) else 1)


if __name__ == "__main__":
    main()
