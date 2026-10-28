#!/usr/bin/env python3
"""
fibre_loop.py - the fibre loop: the fibre loop of OpenAI's identity (2.6) three ways.

The fibre loop is the inner sum of identity (2.6) in "The Exponential Time Hypothesis is
False" (OpenAI, Oct 2026): for a group g of k variables with integer measurement rows C_g,
                B_g(y) = sum_{u in {0,1}^k : C_g u = y}  prod_{e assigned to g} Phi_e(u)
for every attained measurement vector y - the histogram of weighted assignments indexed by
the measurement. The paper's step 2 enumerates all 2^k assignments once per measurement
tuple; the toy in measurement_toy.py enumerates once and buckets by y.

Three methods, all exact modulo p = 2^31 - 1 (a production version would use several primes
and the Chinese remainder theorem; the cost per prime is what is measured), all built on the
SAME vectorised NumPy kernel for enumerating assignments and evaluating clause polynomials,
so that only the algorithm differs:

  bucket    enumerate all 2^k assignments of the group, bucket by y.                [baseline]
  dict      split the group's clause hypergraph into connected components. The measurement
            is additive across components, so the histogram of the group is the convolution
            of the component histograms. Each component is enumerated with the kernel
            (sum_c 2^|c| assignments instead of 2^k) and the convolution is done on sparse
            dictionaries.
  packed    same decomposition; each component histogram is packed into ONE big integer by
            Kronecker substitution (one 80-bit lane per vector of the measurement box), so
            the convolution is one GMP multiplication; lanes are reduced mod p after every
            multiplication (NumPy on the byte image). This is the additive-homomorphic packed
            representation: lanes add under multiplication exactly as measurements add under
            concatenation of partial assignments.

Instances
  --mode design     the group, forms, shifts and measurement rows of measurement_toy's aux
                    design (the paper's Section 2.7 construction at r = 2, Q = 2, R = 2) at
                    t = 1, for a foreign tuple attained by a random assignment.
  --mode synthetic  the within-group clause structure of a random bounded-occurrence 3-CNF
                    under random variable groups and random clause assignments (the paper's
                    Step 4 schedules), with Boolean own inputs, fixed foreign inputs, and a
                    synthetic measurement of c sparse +-1 rows per group. This mimics the
                    paper's regime (few shared rows, c << k), which the real construction only
                    reaches at enormous n.

Every method's histogram is compared with the others and, for small k, with the exact
rational fibre sums of measurement_toy.py. The script prints times, counts (assignments
enumerated, lanes, bits) and peak memory.

Requires numpy and gmpy2; imports measurement_toy.py from the same directory.
"""
from __future__ import annotations

import argparse
import math
import os
import random
import resource
import sys
import time
import tracemalloc
from fractions import Fraction

import numpy as np
import gmpy2
from gmpy2 import mpz

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import measurement_toy as P  # noqa: E402

PRIME = (1 << 31) - 1
LANE_BYTES = 10          # 80-bit lanes: a lane holds < p^2 * 2^|component| after one multiplication


def residue(x) -> int:
    x = Fraction(x)
    return (x.numerator % PRIME) * pow(x.denominator % PRIME, -1, PRIME) % PRIME


# ---------------------------------------------------------------------------
# the vectorised kernel shared by every method
# ---------------------------------------------------------------------------
def histogram_kernel(vs, clauses, rows, chunk_bits=20):
    """Histogram {measurement vector: weight mod p} of the assignments of the variables vs,
    weights = product of the clause polynomials (clauses whose own support lies in vs), with
    the measurement rows restricted to vs. Enumerates 2^|vs| assignments in NumPy chunks."""
    pos = {v: i for i, v in enumerate(vs)}
    vset = set(vs)
    rrows = [{v: c for v, c in row.items() if v in vset} for row in rows]
    dneg = [sum(-c for c in row.values() if c < 0) for row in rrows]
    radix = [dn + sum(c for c in row.values() if c > 0) + 1 for dn, row in zip(dneg, rrows)]
    base = [math.prod(radix[:i]) for i in range(len(rrows))]
    lanes = math.prod(radix)
    k = len(vs)
    hist = np.zeros(lanes, dtype=np.int64)
    total = 1 << k
    step = 1 << min(chunk_bits, k)
    for start in range(0, total, step):
        U = np.arange(start, start + step, dtype=np.int64)
        bits = {}

        def bit(v):
            if v not in bits:
                bits[v] = (U >> pos[v]) & 1
            return bits[v]
        w = np.ones(step, dtype=np.int64)
        for signs, forms, shifts in clauses:
            prod = np.ones(step, dtype=np.int64)
            for s, form, sh in zip(signs, forms, shifts):
                X = np.full(step, sh, dtype=np.int64)
                for v, c in form.items():
                    X += c * bit(v)                    # c < 2^31 and at most k terms: fits int64
                X %= PRIME
                L = X if s else (1 - X) % PRIME
                prod = prod * ((1 - L) % PRIME) % PRIME
            w = w * ((1 - prod) % PRIME) % PRIME
        idx = np.zeros(step, dtype=np.int64)
        for i, row in enumerate(rrows):
            yi = np.zeros(step, dtype=np.int64)
            for v, c in row.items():
                yi += c * bit(v)
            idx += (yi + dneg[i]) * base[i]
        # exact in float64: weights < 2^31, at most 2^20 per lane per chunk
        hist += np.bincount(idx, weights=w.astype(np.float64), minlength=lanes).astype(np.int64)
        hist %= PRIME
    nz = np.nonzero(hist)[0]
    out = {}
    for i in nz:
        j = int(i)
        yv = []
        for rad, dn in zip(radix, dneg):
            yv.append(j % rad - dn)
            j //= rad
        out[tuple(yv)] = int(hist[int(i)])
    return out


# ---------------------------------------------------------------------------
# a group instance: variables, clauses (signs, own forms, shifts), measurement rows
# ---------------------------------------------------------------------------
class GroupInstance:
    def __init__(self, S, clauses, rows, exact_clauses=None, label=""):
        self.S = list(S)
        self.k = len(S)
        self.clauses = clauses                 # (signs, [form dict var->residue], [shift residue])
        self.exact_clauses = exact_clauses     # same with Fractions, for the exact check
        self.rows = rows
        self.label = label
        self.Dneg = [sum(-c for c in row.values() if c < 0) for row in rows]
        self.Dpos = [sum(c for c in row.values() if c > 0) for row in rows]
        self.radix = [a + b + 1 for a, b in zip(self.Dneg, self.Dpos)]
        self.lanes = math.prod(self.radix)
        self.base = [math.prod(self.radix[:i]) for i in range(len(rows))]
        self.m_g = len(clauses)

    @staticmethod
    def from_design(n, r, seed, t=1):
        rng = random.Random(seed)
        cnf = P.random_bounded_occurrence_3cnf(n, n, rng)
        D = P.build_aux_design(cnf, r, 2, 2, rng)
        pt = D.at(t)
        own, lam = D.numeric_forms(t)
        g = max(range(len(pt.groups)), key=lambda h: len(pt.groups[h]))
        S = list(pt.groups[g])
        x = [rng.randrange(2) for _ in range(n)]
        y = [tuple(sum(c * x[v] for v, c in row.items()) for row in pt.rows[h]) for h in range(len(pt.groups))]
        clauses, exact = [], []
        for e in range(len(D.occ)):
            if pt.sigma[e] != g:
                continue
            signs = [s for _, s in D.occ[e]]
            forms, shifts, eforms, eshifts = [], [], [], []
            for kk in range(len(D.occ[e])):
                form = own[e][kk]
                shift = sum((sum((c * yy for c, yy in zip(coefs, y[h])), 0) for h, coefs in lam[e][kk].items()), 0)
                forms.append({v: residue(c) for v, c in form.items()})
                shifts.append(residue(shift))
                eforms.append(dict(form))
                eshifts.append(Fraction(shift))
            clauses.append((signs, forms, shifts))
            exact.append((signs, eforms, eshifts))
        rows = [dict(row) for row in pt.rows[g]]
        return GroupInstance(S, clauses, rows, exact, label=f"aux design r = {r}, t = {t}, group {g + 1}")

    @staticmethod
    def synthetic(n, r, c_rows, seed, support=3):
        rng = random.Random(seed)
        cnf = P.random_bounded_occurrence_3cnf(n, n, rng)
        occ = [[(abs(l) - 1, l > 0) for l in cl] for cl in cnf.clauses]
        # group 0 gets exactly k = n / r random variables; the others are distributed at random
        k = n // r
        S = sorted(rng.sample(range(n), k))
        a = [0] * n
        for v in range(n):
            if v not in set(S):
                a[v] = rng.randrange(1, r)
        d = [rng.randrange(r) for _ in range(len(occ))]
        g = 0
        clauses, exact = [], []
        for e, cl in enumerate(occ):
            if d[e] != g:
                continue
            signs = [s for _, s in cl]
            forms, shifts = [], []
            for v, _ in cl:
                if a[v] == g:
                    forms.append({v: 1})
                    shifts.append(0)
                else:
                    forms.append({})
                    shifts.append(rng.randrange(PRIME))      # a fixed foreign input
            clauses.append((signs, forms, shifts))
            exact.append((signs, [dict(f) for f in forms], [Fraction(s) for s in shifts]))
        rows = []
        for _ in range(c_rows):
            vs = rng.sample(S, min(support, len(S)))
            rows.append({v: rng.choice((-1, 1)) for v in vs})
        return GroupInstance(S, clauses, rows, exact,
                             label=f"synthetic: random schedules at r = {r}, {c_rows} sparse rows of support {support}")

    # ---- decomposition ----
    def components(self):
        parent = {v: v for v in self.S}

        def find(v):
            while parent[v] != v:
                parent[v] = parent[parent[v]]
                v = parent[v]
            return v
        supports = []
        for signs, forms, shifts in self.clauses:
            supp = set()
            for f in forms:
                supp |= set(f)
            supports.append(frozenset(supp))
            vs = list(supp)
            for v in vs[1:]:
                parent[find(v)] = find(vs[0])
        comps = {}
        for v in self.S:
            comps.setdefault(find(v), []).append(v)
        out = []
        const = 1
        assigned = 0
        for vs in comps.values():
            vset = set(vs)
            cs = [c for c, supp in zip(self.clauses, supports) if supp and supp <= vset]
            assigned += len(cs)
            out.append((sorted(vs), cs))
        for (signs, forms, shifts), supp in zip(self.clauses, supports):
            if not supp:                       # all three inputs foreign: a constant factor
                prod = 1
                for s, sh in zip(signs, shifts):
                    L = sh if s else (1 - sh) % PRIME
                    prod = prod * ((1 - L) % PRIME) % PRIME
                const = const * ((1 - prod) % PRIME) % PRIME
                assigned += 1
        assert assigned == len(self.clauses)
        self.const = const
        return sorted(out, key=lambda c: len(c[0]))

    # ---- method 1 ----
    def bucket(self):
        h = histogram_kernel(self.S, self.clauses, self.rows)
        self.stats_bucket = {"enumerated": 1 << self.k, "entries": len(h)}
        return h

    # ---- method 2 ----
    def dict_conv(self):
        comps = self.components()
        acc = {tuple(0 for _ in self.rows): self.const}
        ops = enumerated = 0
        for vs, cs in comps:
            enumerated += 1 << len(vs)
            h = histogram_kernel(vs, cs, self.rows)
            new = {}
            for y1, v1 in acc.items():
                for y2, v2 in h.items():
                    y = tuple(p_ + q_ for p_, q_ in zip(y1, y2))
                    new[y] = (new.get(y, 0) + v1 * v2) % PRIME
                    ops += 1
            acc = new
        self.stats_dict = {"components": len(comps), "largest": max(len(v) for v, _ in comps),
                           "enumerated": enumerated, "dict_ops": ops, "entries": len(acc)}
        return {y: v for y, v in acc.items() if v}

    # ---- method 3 ----
    def _reduce_lanes(self, acc):
        """Reduce every 80-bit lane of the packed integer mod p (NumPy on the byte image)
        and repack. Returns (mpz, lane values as an array)."""
        nbytes = self.lane_bytes
        raw = int(acc).to_bytes(self.lanes * nbytes, "little")
        img = np.frombuffer(raw, dtype=np.uint8).reshape(self.lanes, nbytes).astype(np.int64)
        val = np.zeros(self.lanes, dtype=np.int64)
        for b in range(nbytes - 1, -1, -1):              # Horner over bytes, reduced at every step
            val = (val * 256 + img[:, b]) % PRIME
        out = np.zeros((self.lanes, nbytes), dtype=np.uint8)
        for b in range(4):
            out[:, b] = (val >> (8 * b)) & 255
        return mpz(int.from_bytes(out.tobytes(), "little")), val

    def packed_conv(self):
        comps = self.components()
        # lane width: after a reduction every lane is < p; one multiplication adds at most
        # min(2^|c|, lanes) products of two residues, so < p^2 * min(2^|c|, lanes)
        largest = max(len(v) for v, _ in comps)
        need = 62 + min(largest, max(1, math.ceil(math.log2(self.lanes)))) + 1
        self.lane_bytes = max(LANE_BYTES, (need + 7) // 8)
        W = 8 * self.lane_bytes
        acc = mpz(self.const)
        enumerated = mults = 0
        maxbits = 0
        t_kernel = t_mul = t_reduce = 0.0
        for vs, cs in comps:
            enumerated += 1 << len(vs)
            t0 = time.perf_counter()
            h = histogram_kernel(vs, cs, self.rows)
            t_kernel += time.perf_counter() - t0
            vset = set(vs)
            dneg = [sum(-c for v, c in row.items() if c < 0 and v in vset) for row in self.rows]
            packed = 0
            for yv, v in h.items():
                idx = sum((yi + dn) * b for yi, dn, b in zip(yv, dneg, self.base))
                packed += v << (W * idx)
            t0 = time.perf_counter()
            acc = acc * mpz(packed)
            t_mul += time.perf_counter() - t0
            mults += 1
            maxbits = max(maxbits, int(gmpy2.bit_length(acc)))
            t0 = time.perf_counter()
            acc, val = self._reduce_lanes(acc)
            t_reduce += time.perf_counter() - t0
        _, val = self._reduce_lanes(acc)
        hist = {}
        for i in np.nonzero(val)[0]:
            j = int(i)
            yv = []
            for rad, dn in zip(self.radix, self.Dneg):
                yv.append(j % rad - dn)
                j //= rad
            hist[tuple(yv)] = int(val[int(i)])
        self.stats_packed = {"components": len(comps), "largest": max(len(v) for v, _ in comps),
                             "enumerated": enumerated, "mults": mults, "lane_bits": W, "lanes": self.lanes,
                             "max_bits": maxbits, "t_kernel": t_kernel, "t_mul": t_mul, "t_reduce": t_reduce}
        return hist

    # ---- exact check ----
    def exact_histogram(self):
        hist = {}
        pos = {v: i for i, v in enumerate(self.S)}
        for bits in range(1 << self.k):
            u = {v: (bits >> i) & 1 for v, i in pos.items()}
            w = Fraction(1)
            for signs, forms, shifts in self.exact_clauses:
                inputs = [sh + sum((c * u[v] for v, c in form.items()), 0) for form, sh in zip(forms, shifts)]
                w *= P.clause_poly(signs, inputs, 1)
                if w == 0:
                    break
            if w == 0:
                continue
            yv = tuple(sum(c * u[v] for v, c in row.items()) for row in self.rows)
            hist[yv] = hist.get(yv, 0) + w
        return {y: residue(v) for y, v in hist.items() if residue(v)}


def timed(fn):
    tracemalloc.start()
    t0 = time.perf_counter()
    out = fn()
    dt = time.perf_counter() - t0
    _, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    return out, dt, peak


def run(G: GroupInstance, dict_limit: int = 22, bucket_limit: int = 24):
    k, c = G.k, len(G.rows)
    print(f"{G.label}: k = {k} variables, {G.m_g} clauses assigned, c_g = {c} rows, box = {G.lanes} lanes, "
          f"lanes / 2^k = {G.lanes / (1 << k):.3g}")
    comps = G.components()
    sizes = sorted((len(v) for v, _ in comps), reverse=True)
    print(f"  within-group components: {len(comps)}, sizes {sizes[:10]}{' ...' if len(sizes) > 10 else ''}; "
          f"assignments enumerated by the component methods: {sum(1 << s for s in sizes)} vs 2^k = {1 << k}")
    res = {}
    if k <= bucket_limit:
        h1, t1, m1 = timed(G.bucket)
        res["bucket"] = t1
        print(f"  bucket (2^k enumeration):             {t1*1e3:9.1f} ms   attained y: {len(h1)}")
    else:
        h1 = t1 = None
        print(f"  bucket (2^k enumeration):             skipped (k > {bucket_limit})")
    hd = None
    if k <= dict_limit:
        hd, td, md = timed(G.dict_conv)
        res["dict"] = td
        sd = G.stats_dict
        print(f"  components + dict convolution:         {td*1e3:9.1f} ms   {sd['dict_ops']} dict ops")
    hp, tp, mp = timed(G.packed_conv)
    res["packed"] = tp
    sp = G.stats_packed
    print(f"  components + packed convolution (GMP): {tp*1e3:9.1f} ms   {sp['mults']} multiplications, 80-bit lanes, "
          f"largest product {sp['max_bits']} bits ({sp['max_bits']/8/1e6:.2f} MB); kernel {sp['t_kernel']*1e3:.0f} ms, "
          f"GMP {sp['t_mul']*1e3:.0f} ms, lane reduction {sp['t_reduce']*1e3:.0f} ms; peak Python heap {mp/1e6:.1f} MB")
    ok = (h1 is None or h1 == hp) and (hd is None or hd == hp)
    if k <= 12 and G.exact_clauses is not None:
        he = G.exact_histogram()
        ok = ok and he == hp
        print(f"  check against exact rational fibre sums: {'OK' if he == hp else 'FAIL'}")
    print(f"  histograms agree (mod p): {'OK' if ok else 'FAIL'}"
          + (f";  packed vs bucket: {t1/tp:5.2f}x" if t1 is not None else "")
          + (f";  packed vs dict: {res['dict']/tp:5.2f}x" if hd is not None else ""))
    if not ok:
        raise SystemExit("histograms disagree")
    return res, sizes


def component_statistics(n, seeds=range(3)):
    print(f"\nWithin-group clause hypergraph under random schedules (one point), n = m = {n}:")
    for r in (2, 3, 4, 5, 6, 8, 12):
        rows = []
        for seed in seeds:
            rng = random.Random(seed)
            cnf = P.random_bounded_occurrence_3cnf(n, n, rng)
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
            sizes = {}
            for v in range(n):
                sizes[find(v)] = sizes.get(find(v), 0) + 1
            rows.append(max(sizes.values()))
        print(f"  r = {r:2d}: largest within-group component {min(rows)}..{max(rows)} of a group of {n // r} "
              f"(expected within-group degree 3/r = {3 / r:.2f})")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", choices=["design", "synthetic"], default="synthetic")
    ap.add_argument("--sizes", default="16,18,20,22,24", help="group sizes k (n = r k)")
    ap.add_argument("--r", type=int, default=4)
    ap.add_argument("--rows", type=int, default=8, help="synthetic mode: measurement rows per group")
    ap.add_argument("--support", type=int, default=3, help="synthetic mode: variables per row")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--dict-limit", type=int, default=22, help="largest k for which the dict method is run")
    ap.add_argument("--bucket-limit", type=int, default=24, help="largest k for which the 2^k enumeration is run")
    ap.add_argument("--no-stats", action="store_true")
    args = ap.parse_args()
    print(f"Week 2: the fibre loop of (2.6) three ways on one shared vectorised kernel; arithmetic mod p = 2^31 - 1; "
          f"numpy {np.__version__}, gmpy2 {gmpy2.version()}")
    summary = []
    for k in (int(s) for s in args.sizes.split(",")):
        n = args.r * k
        if args.mode == "design":
            G = GroupInstance.from_design(n, args.r, args.seed)
        else:
            G = GroupInstance.synthetic(n, args.r, args.rows, args.seed, args.support)
        res, sizes = run(G, args.dict_limit, args.bucket_limit)
        summary.append((G.k, len(G.rows), G.lanes, sizes[0], res))
        print()
    print("Summary (ms):   k  rows    lanes  largest comp     bucket       dict     packed   packed/bucket")
    for k, c, lanes, big, res in summary:
        tp = res["packed"] * 1e3
        tb = f"{res['bucket']*1e3:9.1f}" if "bucket" in res else "        -"
        td = f"{res['dict']*1e3:9.1f}" if "dict" in res else "        -"
        ratio = f"{res['bucket']/res['packed']:6.2f}x" if "bucket" in res else "      -"
        print(f"               {k:2d}  {c:4d}  {lanes:7d}  {big:12d}  {tb}  {td}  {tp:9.1f}   {ratio}")
    if not args.no_stats:
        component_statistics(4000)
    print(f"\npeak RSS of this process: {resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e3:.0f} MB")


if __name__ == "__main__":
    main()
