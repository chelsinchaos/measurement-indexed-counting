#!/usr/bin/env python3
"""
constrained_count.py - Week 3 (final) experiment: a component cache keyed by a sparse additive
measurement, for weighted model counting under cardinality / pseudo-Boolean side constraints.

Task. Given a CNF F over n variables, literal weights w0[v], w1[v], and c integer rows
a^(1..c) (a "sparse additive measurement" A x), compute the histogram
        H(y) = sum over models x of F with A x = y  of  prod_v w_{x_v}[v]
for every attained y. With c = 1 and a = (1, ..., 1) this is the weighted count of models by
cardinality; with other integer rows it is a pseudo-Boolean side constraint; with c = 0 it is
plain weighted model counting.

Engine. A miniature sharpSAT-style exact counter: unit propagation, decomposition of the
residual formula into connected components, and a cache keyed by the canonical residual
component. The only change from a standard counter is the VALUE type: a cached component
holds a sparse histogram over the measurement instead of a number, free variables contribute
factors (w0 + w1 z^{a_v}), and components are combined by convolution because the
measurement is additive across disjoint variable sets. This is the same algebra as
evaluating a d-DNNF in the polynomial semiring (generating functions), done top-down with
caching instead of by compilation.

Comparison. The standard route for a cardinality side constraint with an ordinary counter is
to ENCODE it in CNF. The script builds a functional (counting-preserving) sequential counter
for sum x_v = b, appends it to F, and runs the SAME engine with c = 0. Both arms share the
code, the propagation, the decomposition, the cache and the branching heuristic, so the
difference is the treatment of the constraint alone. Every histogram is checked against brute
force for small n, and grid totals against OEIS A006506 (independent sets of the n x n grid).

Standard library only.
"""
from __future__ import annotations

import argparse
import itertools
import math
import random
import sys
import time
from dataclasses import dataclass, field

sys.setrecursionlimit(100000)


# ---------------------------------------------------------------------------
# sparse histograms over Z^c: dict  tuple -> int
# ---------------------------------------------------------------------------
def hconv(a: dict, b: dict) -> dict:
    if len(a) == 1 and next(iter(a)) == ():
        pass
    out: dict = {}
    for y1, v1 in a.items():
        for y2, v2 in b.items():
            y = tuple(p + q for p, q in zip(y1, y2))
            out[y] = out.get(y, 0) + v1 * v2
    return out


def hadd(a: dict, b: dict) -> dict:
    out = dict(a)
    for y, v in b.items():
        out[y] = out.get(y, 0) + v
    return out


# ---------------------------------------------------------------------------
# the counter
# ---------------------------------------------------------------------------
@dataclass
class Stats:
    nodes: int = 0
    cache_hits: int = 0
    cache_misses: int = 0
    conflicts: int = 0
    max_component: int = 0
    hist_entries_cached: int = 0
    convolutions: int = 0


class Counter:
    """Exact weighted counting with a histogram-valued component cache."""

    def __init__(self, n: int, clauses, w0=None, w1=None, rows=None, order: str = "dynamic"):
        self.n = n
        self.clauses = tuple(tuple(sorted(set(c))) for c in clauses)
        self.w0 = list(w0) if w0 is not None else [1] * n
        self.w1 = list(w1) if w1 is not None else [1] * n
        rows = rows or []
        self.c = len(rows)
        self.vec = [tuple(int(row[v]) for row in rows) for v in range(n)]
        self.zero = tuple(0 for _ in range(self.c))
        self.order = order
        self.cache: dict = {}
        self.stats = Stats()

    # literal / free-variable factors
    def lit_factor(self, lit: int) -> dict:
        v = abs(lit) - 1
        if lit > 0:
            return {self.vec[v]: self.w1[v]}
        return {self.zero: self.w0[v]}

    def var_factor(self, v: int) -> dict:
        h = {self.zero: self.w0[v]}
        h[self.vec[v]] = h.get(self.vec[v], 0) + self.w1[v]
        return h

    # unit propagation on a clause list; returns (forced literals, residual clauses) or None on conflict
    @staticmethod
    def propagate(clauses, lits):
        forced = []
        assign = {}
        queue = list(lits)
        cur = list(clauses)
        while queue:
            lit = queue.pop()
            v = abs(lit)
            if v in assign:
                if assign[v] != (lit > 0):
                    return None
                continue
            assign[v] = lit > 0
            forced.append(lit)
            nxt = []
            for cl in cur:
                if lit in cl:
                    continue
                if -lit in cl:
                    cl = tuple(x for x in cl if x != -lit)
                    if not cl:
                        return None
                    if len(cl) == 1:
                        queue.append(cl[0])
                nxt.append(cl)
            cur = nxt
        # a clause may have become unit while its literal is already queued/assigned: clean up
        res = []
        for cl in cur:
            if any(abs(x) in assign for x in cl):
                # cannot happen: assigned literals were removed above
                raise AssertionError
            res.append(cl)
        return forced, tuple(res)

    @staticmethod
    def components(clauses):
        parent = {}

        def find(x):
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x
        for cl in clauses:
            vs = [abs(x) for x in cl]
            for v in vs:
                parent.setdefault(v, v)
            for v in vs[1:]:
                parent[find(v)] = find(vs[0])
        groups: dict = {}
        for cl in clauses:
            groups.setdefault(find(abs(cl[0])), []).append(cl)
        return list(groups.values())

    def choose(self, clauses):
        if self.order == "static":
            return min(abs(x) for cl in clauses for x in cl)
        occ: dict = {}
        for cl in clauses:
            for x in cl:
                occ[abs(x)] = occ.get(abs(x), 0) + 1
        return max(occ, key=lambda v: (occ[v], -v))

    def solve(self, clauses) -> dict:
        """Histogram of the models of a connected clause set over exactly its variables."""
        key = tuple(sorted(clauses))
        hit = self.cache.get(key)
        if hit is not None:
            self.stats.cache_hits += 1
            return hit
        self.stats.cache_misses += 1
        self.stats.nodes += 1
        vars_here = {abs(x) for cl in clauses for x in cl}
        self.stats.max_component = max(self.stats.max_component, len(vars_here))
        v = self.choose(clauses)
        total: dict = {}
        for lit in (v, -v):
            prop = self.propagate(clauses, [lit])
            if prop is None:
                self.stats.conflicts += 1
                continue
            forced, residual = prop
            h = {self.zero: 1}
            for f in forced:
                h = hconv(h, self.lit_factor(f))
            left = vars_here - {abs(f) for f in forced}
            if residual:
                comps = self.components(residual)
                in_res = {abs(x) for cl in residual for x in cl}
                for fv in left - in_res:
                    h = hconv(h, self.var_factor(fv - 1))
                for comp in comps:
                    h = hconv(h, self.solve(tuple(comp)))
                    self.stats.convolutions += 1
            else:
                for fv in left:
                    h = hconv(h, self.var_factor(fv - 1))
            total = hadd(total, h)
        self.cache[key] = total
        self.stats.hist_entries_cached += len(total)
        return total

    def run(self) -> dict:
        """Histogram over all n variables."""
        units = [cl[0] for cl in self.clauses if len(cl) == 1]
        prop = self.propagate(self.clauses, units) if units else ([], self.clauses)
        if prop is None:
            return {}
        forced, residual = prop
        h = {self.zero: 1}
        for f in forced:
            h = hconv(h, self.lit_factor(f))
        in_res = {abs(x) for cl in residual for x in cl}
        fset = {abs(f) for f in forced}
        for v in range(1, self.n + 1):
            if v not in in_res and v not in fset:
                h = hconv(h, self.var_factor(v - 1))
        for comp in self.components(residual):
            h = hconv(h, self.solve(tuple(comp)))
        return {y: w for y, w in h.items() if w}


# ---------------------------------------------------------------------------
# the standard route: a functional sequential counter for  sum_{v in V} x_v = b
# ---------------------------------------------------------------------------
def encode_cardinality(n: int, V, b: int, clauses):
    """Append clauses defining s(i, j) <-> (at least j of x_{V[0..i-1]} are true) for
    1 <= j <= min(i, b + 1), and assert s(m, b) and not s(m, b + 1). Every auxiliary is a
    function of x, so the model count of the extended formula equals the number of models
    of F with sum x_V = b. Returns (new n, new clause list, number of auxiliaries)."""
    V = list(V)
    m = len(V)
    out = list(clauses)
    nv = n
    if b > m:
        out.append((1,)); out.append((-1,))                  # unsatisfiable
        return nv, out, 0
    J = min(b + 1, m)
    s = {}                                                  # (i, j) -> variable; absent = FALSE (j > i)
    for i in range(1, m + 1):
        x = V[i - 1]
        for j in range(1, min(i, J) + 1):
            nv += 1
            p = nv
            s[(i, j)] = p
            q = s.get((i - 1, j))                           # FALSE when j > i - 1
            if j == 1:                                      # r = s(i-1, 0) = TRUE:  p <-> q or x
                if q is None:                               # p <-> x
                    out.append((-p, x)); out.append((-x, p))
                else:
                    out.append((-p, q, x)); out.append((-q, p)); out.append((-x, p))
            else:
                r = s.get((i - 1, j - 1))                   # variable (j - 1 <= i - 1 always holds here)
                assert r is not None
                if q is None:                               # p <-> x and r
                    out.append((-p, x)); out.append((-p, r)); out.append((-x, -r, p))
                else:                                       # p <-> q or (x and r)
                    out.append((-p, q, x)); out.append((-p, q, r)); out.append((-q, p)); out.append((-x, -r, p))
    if b >= 1:
        out.append((s[(m, b)],))
    if b + 1 <= m:
        out.append((-s[(m, b + 1)],))
    return nv, out, nv - n


# ---------------------------------------------------------------------------
# instances
# ---------------------------------------------------------------------------
def grid_independent_sets(W: int, H: int):
    """Variables = vertices (row-major), clauses = one 2-clause per edge."""
    def vid(r, c):
        return r * W + c + 1
    clauses = []
    for r in range(H):
        for c in range(W):
            if c + 1 < W:
                clauses.append((-vid(r, c), -vid(r, c + 1)))
            if r + 1 < H:
                clauses.append((-vid(r, c), -vid(r + 1, c)))
    return W * H, clauses


def random_3cnf(n: int, m: int, rng: random.Random):
    return [tuple(v if rng.random() < 0.5 else -v for v in rng.sample(range(1, n + 1), 3)) for _ in range(m)]


def brute_force(n, clauses, w0, w1, rows):
    hist: dict = {}
    for bits in range(1 << n):
        x = [(bits >> i) & 1 for i in range(n)]
        if all(any((x[abs(l) - 1] == 1) == (l > 0) for l in cl) for cl in clauses):
            y = tuple(sum(row[v] * x[v] for v in range(n)) for row in rows)
            w = 1
            for v in range(n):
                w *= w1[v] if x[v] else w0[v]
            hist[y] = hist.get(y, 0) + w
    return hist


OEIS_A006506 = {1: 2, 2: 7, 3: 63, 4: 1234, 5: 55447, 6: 5598861, 7: 1280128950, 8: 660647962955}


def fmt_stats(st: Stats, dt: float) -> str:
    return (f"{dt*1e3:9.1f} ms, {st.nodes} nodes, cache {st.cache_misses} entries / {st.cache_hits} hits, "
            f"largest component {st.max_component} vars, {st.hist_entries_cached} histogram entries cached")


def experiment_grid(sizes, order, encode_sizes):
    print("== Grid independent sets by size: histogram cache vs cardinality encoding ==")
    for W, H in sizes:
        n, clauses = grid_independent_sets(W, H)
        ones = [[1] * n]
        C = Counter(n, clauses, rows=ones, order=order)
        t0 = time.perf_counter()
        hist = C.run()
        dt = time.perf_counter() - t0
        total = sum(hist.values())
        mode = max(hist, key=lambda y: hist[y])[0]
        check = ""
        if W == H and W in OEIS_A006506:
            check = f"; total {total} {'= OEIS A006506' if total == OEIS_A006506[W] else 'DIFFERS from OEIS A006506 -> FAIL'}"
        if n <= 20:
            bf = brute_force(n, clauses, C.w0, C.w1, ones)
            check += f"; brute force {'OK' if bf == hist else 'FAIL'}"
        print(f"{W}x{H} grid, n = {n}: histogram of {len(hist)} sizes in one run: {fmt_stats(C.stats, dt)}{check}")
        print(f"      sizes 0..{max(y[0] for y in hist)}: {[hist.get((k,), 0) for k in range(min(len(hist), 12))]}{' ...' if len(hist) > 12 else ''}; mode b = {mode}")
        if (W, H) in encode_sizes:
            for b in sorted({mode, mode // 2, 1}):
                n2, cl2, aux = encode_cardinality(n, range(1, n + 1), b, clauses)
                C2 = Counter(n2, cl2, order=order)
                t0 = time.perf_counter()
                h2 = C2.run()
                dt2 = time.perf_counter() - t0
                val = h2.get((), 0)
                ok = val == hist.get((b,), 0)
                print(f"      encoding, b = {b:2d}: {aux} auxiliaries, {len(cl2)} clauses: {fmt_stats(C2.stats, dt2)}; "
                      f"count {val} {'OK' if ok else 'FAIL'}; this one b cost {dt2/dt:5.2f}x the whole histogram")
        print()


def experiment_random(n, m, seed, order, bs, weighted):
    rng = random.Random(seed)
    clauses = random_3cnf(n, m, rng)
    w0 = [1] * n
    w1 = [rng.randint(1, 3) for _ in range(n)] if weighted else [1] * n
    ones = [[1] * n]
    print(f"== Random 3-CNF n = {n}, m = {m}{', literal weights w1 in 1..3' if weighted else ''}: histogram cache vs encoding ==")
    C = Counter(n, clauses, w0, w1, rows=ones, order=order)
    t0 = time.perf_counter()
    hist = C.run()
    dt = time.perf_counter() - t0
    total = sum(hist.values())
    print(f"histogram of {len(hist)} cardinalities in one run: {fmt_stats(C.stats, dt)}; weighted total {total}")
    if n <= 20:
        bf = brute_force(n, clauses, w0, w1, ones)
        print(f"      brute force: {'OK' if bf == hist else 'FAIL'}")
    C0 = Counter(n, clauses, w0, w1, order=order)
    t0 = time.perf_counter()
    plain = C0.run()
    dt0 = time.perf_counter() - t0
    print(f"      plain weighted count (no measurement, same engine): {fmt_stats(C0.stats, dt0)}; "
          f"total {plain.get((), 0)} {'OK' if plain.get((), 0) == total else 'FAIL'}; histogram cost {dt/dt0:5.2f}x the plain count")
    for b in bs:
        n2, cl2, aux = encode_cardinality(n, range(1, n + 1), b, clauses)
        C2 = Counter(n2, cl2, w0 + [1] * (n2 - n), w1 + [1] * (n2 - n), order=order)
        t0 = time.perf_counter()
        h2 = C2.run()
        dt2 = time.perf_counter() - t0
        val = h2.get((), 0)
        ok = val == hist.get((b,), 0)
        print(f"      encoding, b = {b:2d}: {aux} auxiliaries, {len(cl2)} clauses: {fmt_stats(C2.stats, dt2)}; "
              f"count {val} {'OK' if ok else 'FAIL'}; this one b cost {dt2/dt:5.2f}x the whole histogram")
    print()


def experiment_multirow(W, H, order, seed):
    rng = random.Random(seed)
    n, clauses = grid_independent_sets(W, H)
    print(f"== {W}x{H} grid: pseudo-Boolean and two-row measurements ==")
    pb = [[rng.randint(1, 3) for _ in range(n)]]
    C = Counter(n, clauses, rows=pb, order=order)
    t0 = time.perf_counter(); hist = C.run(); dt = time.perf_counter() - t0
    check = ""
    if n <= 20:
        check = f"; brute force {'OK' if brute_force(n, clauses, C.w0, C.w1, pb) == hist else 'FAIL'}"
    print(f"pseudo-Boolean row a_v in 1..3: histogram of {len(hist)} values 0..{max(y[0] for y in hist)}: {fmt_stats(C.stats, dt)}{check}")
    left = [1 if (v % W) < W // 2 else 0 for v in range(n)]
    right = [1 - l for l in left]
    two = [left, right]
    C = Counter(n, clauses, rows=two, order=order)
    t0 = time.perf_counter(); hist = C.run(); dt = time.perf_counter() - t0
    check = ""
    if n <= 20:
        check = f"; brute force {'OK' if brute_force(n, clauses, C.w0, C.w1, two) == hist else 'FAIL'}"
    print(f"two rows (left-half count, right-half count): histogram of {len(hist)} pairs: {fmt_stats(C.stats, dt)}{check}")
    print()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--order", choices=["static", "dynamic"], default="static")
    ap.add_argument("--grids", default="4x4,5x5,6x6,7x7,8x8")
    ap.add_argument("--encode-grids", default="4x4,5x5,6x6")
    ap.add_argument("--random-n", type=int, default=40)
    ap.add_argument("--random-m", type=int, default=100)
    ap.add_argument("--seed", type=int, default=1)
    args = ap.parse_args()
    grids = [tuple(int(t) for t in s.split("x")) for s in args.grids.split(",")]
    enc = {tuple(int(t) for t in s.split("x")) for s in args.encode_grids.split(",")}
    print(f"constrained_count.py: histogram-valued component cache, branching order = {args.order}\n")
    experiment_grid(grids, args.order, enc)
    experiment_multirow(6, 6, args.order, args.seed)
    experiment_random(20, 50, args.seed, "dynamic", [5, 10], weighted=False)
    experiment_random(args.random_n, args.random_m, args.seed, "dynamic", [args.random_n // 2, args.random_n // 4], weighted=True)


if __name__ == "__main__":
    main()
