#!/usr/bin/env python3
"""
measurement_toy.py - an exact, executable model of the #3-SAT counting algorithm of
"The Exponential Time Hypothesis is False" (OpenAI, October 2026), Theorem 1.1, at toy scale.

Implements the clause polynomials, replacement linear forms, measurement designs
(Definition 2.1 of the paper), the measurement-tuple identity (2.6), the restricted sums and
inclusion-exclusion identity (2.9), Hermite interpolation of H from Taylor data, and three
hand-built designs including the paper's own Section 2.7 construction at small parameters.
Every stage is checked against brute-force enumeration. Exact arithmetic throughout
(Python integers and Fractions). Standard library only, Python 3.10+.

This file is the toy extracted from pnp_reference.py (McElveen, 2026) for the paper
"Measurement-Indexed Counting"; the two are identical in the functions they share.

    python3 measurement_toy.py --n 12 --m 12 --q 1              # identities (2.6), (2.9), H(0) = #F
    python3 measurement_toy.py --n 12 --m 12 --design aux       # the paper's Section 2.7 design at r=2, Q=2, R=2
"""
from __future__ import annotations

import argparse
import math
import random
import sys
import time
from dataclasses import dataclass, field
from fractions import Fraction
from typing import Optional, Sequence


@dataclass
class CNF:
    n: int
    clauses: list[tuple[int, ...]]
    pos: list[int] = field(default_factory=list)   # per-clause mask of positive literals
    neg: list[int] = field(default_factory=list)   # per-clause mask of negative literals

    def __post_init__(self):
        self.clauses = [tuple(c) for c in self.clauses]
        for c in self.clauses:
            if any(abs(l) < 1 or abs(l) > self.n or l == 0 for l in c):
                raise ValueError(f"literal out of range in clause {c}")
        self.pos, self.neg = [], []
        for c in self.clauses:
            p = q = 0
            for l in c:
                if l > 0:
                    p |= 1 << (l - 1)
                else:
                    q |= 1 << (-l - 1)
            self.pos.append(p)
            self.neg.append(q)

    @property
    def m(self) -> int:
        return len(self.clauses)

    def satisfied_by(self, B: int) -> bool:
        """Evaluate at the assignment with index B (variable i <-> bit i-1)."""
        notB = ~B
        for p, q in zip(self.pos, self.neg):
            if not (B & p) and not (notB & q):
                return False
        return True

    @staticmethod
    def from_dimacs(text: str) -> "CNF":
        n = 0
        clauses: list[list[int]] = []
        cur: list[int] = []
        for line in text.splitlines():
            line = line.strip()
            if not line or line[0] in "c%":
                continue
            if line.startswith("p"):
                parts = line.split()
                n = int(parts[2])
                continue
            for tok in line.split():
                v = int(tok)
                if v == 0:
                    if cur:
                        clauses.append(cur)
                    cur = []
                else:
                    cur.append(v)
        if cur:
            clauses.append(cur)
        if n == 0:
            n = max(abs(l) for c in clauses for l in c)
        return CNF(n, clauses)

    @staticmethod
    def parse_inline(n: int, text: str) -> "CNF":
        """'1 2 -3 ; -1 4 ; 2' -> three clauses."""
        clauses = [[int(t) for t in part.split()] for part in text.split(";") if part.strip()]
        return CNF(n, clauses)

    @staticmethod
    def random_k_cnf(n: int, m: int, k: int, rng: random.Random, allow_repeats: bool = False) -> "CNF":
        clauses = []
        for _ in range(m):
            if allow_repeats:
                vs = [rng.randint(1, n) for _ in range(k)]
            else:
                vs = rng.sample(range(1, n + 1), k)
            clauses.append([v if rng.random() < 0.5 else -v for v in vs])
        return CNF(n, clauses)

    @staticmethod
    def planted(n: int, m: int, k: int, rng: random.Random) -> tuple["CNF", int]:
        """Random k-CNF satisfied by a hidden random assignment (returned as an index)."""
        hidden = rng.getrandbits(n)
        clauses = []
        while len(clauses) < m:
            vs = rng.sample(range(1, n + 1), k)
            cl = [v if rng.random() < 0.5 else -v for v in vs]
            if any(((hidden >> (abs(l) - 1)) & 1) == (1 if l > 0 else 0) for l in cl):
                clauses.append(cl)
        return CNF(n, clauses), hidden

    @staticmethod
    def forced(n: int, B: int) -> "CNF":
        """F_b: one unit clause per variable; the unique model is Dec(B)."""
        return CNF(n, [[(i + 1) if (B >> i) & 1 else -(i + 1)] for i in range(n)])


def count_solutions(cnf: CNF) -> int:
    """Exact count by enumeration; intended for n <= 22."""
    return sum(1 for B in range(1 << cnf.n) if cnf.satisfied_by(B))


# =============================================================================
# 12. Toy of OpenAI Theorem 1.1: identities (2.6) and (2.9) with hand-built designs
# =============================================================================
#
# Notation follows the paper (Section 2). A 3-CNF in which every variable occurs in at
# most three clauses has clause polynomials Phi_e(X) = 1 - prod_v (1 - L_{e,v}(X_v)) and
# #F = sum_x prod_e Phi_e(x). Replacement linear forms P_{e,v}(T) with P(0) = e_v define
# H(T) = sum_x prod_e Phi_e(<P_{e,v}(T), x>), so H(0) = #F and deg H <= M degP, where M is
# the number of literal occurrences.
#
# At a point t where an (r, tau, D)-design exists (Definition 2.1: groups S_{t,g}, clause
# assignments sigma_e(t), integer measurement matrices C_{t,g}), identity (2.6) computes H(t)
# as a sum over measurement tuples y of prod_g B_g(y), each B_g a sum over the fibre
# {u in {0,1}^{S_g} : C_g u = y_g}. Identity (2.9) obtains the Taylor coefficients of H at t
# from restricted sums Z_{t,E,U} by inclusion-exclusion, and Hermite interpolation over the
# points t = 1..Q recovers H(0).
#
# Three hand-built designs:
#   direct   : C_{t,h} reads every boundary variable of group h directly (one unit row per
#              variable of S_h occurring in a clause assigned elsewhere); groups and rows do
#              not depend on t, so the design holds at every t. With the own-group deformation
#              P = e_v + T w (degP = 1), H is a nonconstant polynomial and the full (2.9) +
#              interpolation pipeline runs with q < n at n = 10..14.
#   collapse : one all-ones row per group; P = e_v (1 - T) + s_h T collapses every foreign
#              occurrence to the group weight at t = 1. Q = 1 forces q = M (exhaustive).
#   aux      : the paper's own construction (Section 2.7) at small parameters: auxiliary
#              coordinates B = ([r] x [r])^Q, swap matrices Pi_s, M_0 = V_0 f(W_0) with f
#              interpolating (-1)^j on an interval of R consecutive integers, measurement rows
#              from the unfixed eigenvectors, rows P_u(T) from (2.16), schedules a_v, d_e and
#              copies from Step 4, and group-split measurement matrices C_{t,g}. The design
#              is certified at every t in [Q]; groups, assignments and rows all depend on t.
#              Rows are genuine +-1 sums of several variables shared by all clauses, which is
#              the mechanism the paper's saving rests on.
#
# Everything is exact (Python integers and Fractions). Every stage is checked against brute
# force. The toy demonstrates correctness and cost, not savings: the saving regime needs
# R << Q with D = r^{2Q} << n (its Section 2.6), which exists only for astronomically large n.

class Poly:
    """Truncated power series in z modulo z^J, coefficients Python ints or Fractions."""
    __slots__ = ("c", "J")

    def __init__(self, c, J: int):
        c = list(c)[:J]
        while c and c[-1] == 0:
            c.pop()
        self.c = c
        self.J = J

    @staticmethod
    def const(a, J: int) -> "Poly":
        return Poly([a], J)

    def is_zero(self) -> bool:
        return not self.c

    def __add__(self, o):
        if not isinstance(o, Poly):
            o = Poly.const(o, self.J)
        n = max(len(self.c), len(o.c))
        return Poly([(self.c[i] if i < len(self.c) else 0) + (o.c[i] if i < len(o.c) else 0) for i in range(n)], self.J)

    __radd__ = __add__

    def __neg__(self):
        return Poly([-a for a in self.c], self.J)

    def __sub__(self, o):
        return self + (-o)

    def __rsub__(self, o):
        return (-self) + o

    def __mul__(self, o):
        if not isinstance(o, Poly):
            return Poly([a * o for a in self.c], self.J)
        if not self.c or not o.c:
            return Poly([], self.J)
        out = [0] * min(self.J, len(self.c) + len(o.c) - 1)
        for i, a in enumerate(self.c):
            if a == 0:
                continue
            for j, b in enumerate(o.c):
                if i + j >= self.J:
                    break
                out[i + j] += a * b
        return Poly(out, self.J)

    __rmul__ = __mul__

    def coeffs(self) -> list:
        return self.c + [0] * (self.J - len(self.c))


def clause_poly(signs: Sequence[bool], inputs: Sequence, one):
    """Phi_e = 1 - prod_v (1 - L_v), L_v = X_v (positive literal) or 1 - X_v (negative)."""
    prod = one
    for s, X in zip(signs, inputs):
        L = X if s else (one - X)
        prod = prod * (one - L)
    return one - prod


def random_bounded_occurrence_3cnf(n: int, m: int, rng: random.Random) -> CNF:
    """Random 3-CNF with distinct variables per clause and every variable in <= 3 clauses."""
    counts = [0] * n
    clauses = []
    for _ in range(m):
        avail = [v for v in range(n) if counts[v] < 3]
        if len(avail) < 3:
            break
        vs = rng.sample(avail, 3)
        for v in vs:
            counts[v] += 1
        clauses.append([(v + 1) if rng.random() < 0.5 else -(v + 1) for v in vs])
    return CNF(n, clauses)


def _poly_eval(coefs: Sequence, t: int):
    return sum((c * t ** j for j, c in enumerate(coefs)), 0)


def _pmul(a: Sequence, b: Sequence) -> list:
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] += x * y
    return out


def _padd(a: Sequence, b: Sequence) -> list:
    n = max(len(a), len(b))
    return [(a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0) for i in range(n)]


def _taylor_shift(coefs: Sequence, t: int) -> list:
    """Coefficients in z of p(t + z) for p given by coefs in T."""
    out = [0] * len(coefs)
    for j, c in enumerate(coefs):
        if c == 0:
            continue
        for i in range(j + 1):
            out[i] += c * math.comb(j, i) * t ** (j - i)
    return out


def _lagrange_coeffs(xs: Sequence[int], ys: Sequence[int]) -> list:
    """Coefficients (low to high, Fractions) of the polynomial through the points (xs, ys)."""
    coefs = [Fraction(0)] * len(xs)
    for i, xi in enumerate(xs):
        num = [Fraction(1)]
        den = Fraction(1)
        for j, xj in enumerate(xs):
            if j != i:
                num = _pmul(num, [Fraction(-xj), Fraction(1)])
                den *= xi - xj
        for k, c in enumerate(num):
            coefs[k] += Fraction(ys[i]) * c / den
    return coefs


def _solve_in_row_space(vec: dict, rows: list) -> Optional[list]:
    """Coefficients mu with sum_i mu_i rows[i] == vec (dicts key -> number), or None if vec lies
    outside the row space. Exact Gaussian elimination over Q on a small system."""
    if not rows:
        return [] if all(c == 0 for c in vec.values()) else None
    keys = sorted(set(vec) | {k for row in rows for k in row})
    A = [[Fraction(row.get(k, 0)) for row in rows] + [Fraction(vec.get(k, 0))] for k in keys]
    ncol = len(rows)
    rix = 0
    pivots = []
    for col in range(ncol):
        p = next((i for i in range(rix, len(A)) if A[i][col] != 0), None)
        if p is None:
            continue
        A[rix], A[p] = A[p], A[rix]
        pv = A[rix][col]
        A[rix] = [x / pv for x in A[rix]]
        for i in range(len(A)):
            if i != rix and A[i][col] != 0:
                f = A[i][col]
                A[i] = [a - f * b for a, b in zip(A[i], A[rix])]
        pivots.append(col)
        rix += 1
    for i in range(rix, len(A)):
        if A[i][-1] != 0:
            return None
    mu = [Fraction(0)] * ncol
    for i, col in enumerate(pivots):
        mu[col] = A[i][-1]
    return mu


@dataclass
class DesignPoint:
    """An (r, tau, D)-design at one point t: groups S_{t,g}, clause assignments sigma_e(t) and
    integer measurement rows of C_{t,g} (each row a dict variable -> coefficient)."""
    groups: list
    group_of: list
    sigma: list
    rows: list


@dataclass
class ToyDesign:
    kind: str
    r: int
    Q: int                      # points at which the design is certified: t = 1..Q
    degP: int                   # degree bound of the replacement forms
    occ: list                   # clause -> [(variable, positive?)]
    P: list                     # P[e][k]: variable -> coefficient list in T (ints or Fractions)
    points: dict                # t -> DesignPoint; key None = the same design at every t
    M: int
    info: dict = field(default_factory=dict)

    def at(self, t: int) -> DesignPoint:
        if t in self.points:
            return self.points[t]
        if None in self.points:
            return self.points[None]
        raise ValueError(f"no design is certified at t = {t}")

    def numeric_forms(self, t: int):
        """Forms at the point t split into the own-group part and, for every other group h, the
        coefficients lambda_{e,v,h} over the rows of C_{t,h} (condition (4), eq. (2.5)).
        Raises if some form's restriction to a foreign group is outside the row space."""
        pt = self.at(t)
        own, lam = [], []
        for e, occs in enumerate(self.occ):
            g = pt.sigma[e]
            own_e, lam_e = [], []
            for k, (v, _s) in enumerate(occs):
                vec = {i: _poly_eval(cs, t) for i, cs in self.P[e][k].items()}
                own_e.append({i: c for i, c in vec.items() if pt.group_of[i] == g and c != 0})
                lam_ek = {}
                for h in range(len(pt.groups)):
                    if h == g:
                        continue
                    rest = {i: c for i, c in vec.items() if pt.group_of[i] == h and c != 0}
                    if not rest:
                        continue
                    mu = _solve_in_row_space(rest, pt.rows[h])
                    if mu is None:
                        raise ValueError(f"design condition (4) fails at t = {t}: clause {e + 1}, "
                                         f"variable {v + 1}, group {h + 1}")
                    lam_ek[h] = mu
                lam_e.append(lam_ek)
            own.append(own_e)
            lam.append(lam_e)
        return own, lam


def _occurrences(cnf: CNF) -> list:
    occ = [[(abs(l) - 1, l > 0) for l in c] for c in cnf.clauses]
    for v in range(cnf.n):
        if sum(1 for c in occ for (w, _) in c if w == v) > 3:
            raise ValueError("variable occurs in more than three clauses; the toy needs the bounded-occurrence form")
    return occ


def build_direct_design(cnf: CNF, r: int, deform: bool, rng: random.Random) -> ToyDesign:
    n = cnf.n
    occ = _occurrences(cnf)
    size = -(-n // r)
    groups = [S for S in (list(range(g * size, min(n, (g + 1) * size))) for g in range(r)) if S]
    r = len(groups)
    group_of = [0] * n
    for g, S in enumerate(groups):
        for v in S:
            group_of[v] = g
    sigma = []
    for c in occ:
        cnt = [0] * r
        for v, _ in c:
            cnt[group_of[v]] += 1
        sigma.append(max(range(r), key=lambda g: (cnt[g], -g)))
    boundary = [set() for _ in range(r)]
    for e, c in enumerate(occ):
        for v, _ in c:
            if group_of[v] != sigma[e]:
                boundary[group_of[v]].add(v)
    rows = [[{v: 1} for v in sorted(boundary[h])] for h in range(r)]
    P = []
    for e, c in enumerate(occ):
        Pe = []
        for v, _ in c:
            form = {v: [1]}
            if deform and group_of[v] == sigma[e]:
                for w in groups[sigma[e]]:
                    coef = rng.choice((-1, 0, 1))
                    if coef:
                        form[w] = [1 if w == v else 0, coef]
            Pe.append(form)
        P.append(Pe)
    M = sum(len(c) for c in occ)
    pt = DesignPoint(groups, group_of, sigma, rows)
    return ToyDesign("direct", r, 0, 1 if deform else 0, occ, P, {None: pt}, M,
                     {"boundary": sum(len(b) for b in boundary), "deform": deform})


def build_collapse_design(cnf: CNF, r: int) -> ToyDesign:
    n = cnf.n
    occ = _occurrences(cnf)
    size = -(-n // r)
    groups = [S for S in (list(range(g * size, min(n, (g + 1) * size))) for g in range(r)) if S]
    r = len(groups)
    group_of = [0] * n
    for g, S in enumerate(groups):
        for v in S:
            group_of[v] = g
    sigma = []
    for c in occ:
        cnt = [0] * r
        for v, _ in c:
            cnt[group_of[v]] += 1
        sigma.append(max(range(r), key=lambda g: (cnt[g], -g)))
    rows = [[{v: 1 for v in groups[h]}] for h in range(r)]
    P = []
    for e, c in enumerate(occ):
        Pe = []
        for v, _ in c:
            if group_of[v] == sigma[e]:
                Pe.append({v: [1]})
            else:
                form = {w: [0, 1] for w in groups[group_of[v]]}       # s_h T
                form[v] = [1, 0]                                      # + e_v (1 - T)
                Pe.append(form)
        P.append(Pe)
    M = sum(len(c) for c in occ)
    return ToyDesign("collapse", r, 1, 1, occ, P, {None: DesignPoint(groups, group_of, sigma, rows)}, M, {})


def build_aux_design(cnf: CNF, r: int, Q: int, R: int, rng: random.Random) -> ToyDesign:
    """Section 2.7 of the paper at small (r, Q, R). Every step is checked as it is built."""
    import itertools
    n = cnf.n
    occ = _occurrences(cnf)
    M = sum(len(c) for c in occ)
    if r < 2 or Q < 1 or not (1 <= R <= Q + 1):
        raise ValueError("aux design needs r >= 2, Q >= 1 and 1 <= R <= Q + 1")

    # Step 1: auxiliary coordinates B = ([r] x [r])^Q, index u = ((a_u(1), d_u(1)), ..., (a_u(Q), d_u(Q))).
    pairs = [(a, d) for a in range(r) for d in range(r)]
    B = list(itertools.product(pairs, repeat=Q))
    Dn = len(B)

    def swap(u, s):
        return u[:s] + ((u[s][1], u[s][0]),) + u[s + 1:]

    def swap_all(u):
        return tuple((d, a) for a, d in u)

    # Step 2: f of degree <= R - 1 with f(j) = (-1)^j on I = {A0, ..., A0 + R - 1}, then M0 = V0 f(W0).
    p = Fraction(r - 1, 2 * r)
    A0 = int(math.floor(Q * p - Fraction(R - 1, 2)))
    A0 = max(0, min(Q - R + 1, A0))
    I = list(range(A0, A0 + R))
    f_coefs = _lagrange_coeffs(I, [(-1) ** j for j in I])

    def apply_W0(vec):
        out = {}
        for u, c in vec.items():
            for s in range(Q):
                w = swap(u, s)
                out[u] = out.get(u, 0) + c / 2
                out[w] = out.get(w, 0) - c / 2
        return {k: v for k, v in out.items() if v != 0}

    def apply_M0(vec):                        # M0 = V0 f(W0), f(W0) by Horner's rule
        res = {}
        for coef in reversed(f_coefs):
            res = apply_W0(res)
            for u, c in vec.items():
                res[u] = res.get(u, 0) + coef * c
            res = {k: v for k, v in res.items() if v != 0}
        return {swap_all(u): c for u, c in res.items()}

    M0_cols = {u: apply_M0({u: Fraction(1)}) for u in B}
    M0_rows = {u: {} for u in B}
    for w, col in M0_cols.items():
        for u, c in col.items():
            M0_rows[u][w] = c

    # Eigenbasis of the swaps: per position e_aa (+1), e_ab + e_ba (+1), e_ab - e_ba (-1).
    per_pos = [({(a, a): 1}, 1) for a in range(r)]
    for a in range(r):
        for b in range(a + 1, r):
            per_pos.append(({(a, b): 1, (b, a): 1}, 1))
            per_pos.append(({(a, b): 1, (b, a): -1}, -1))
    C_rows = []                     # measurement rows: transposes of the eigenvectors with j not in I
    fixed = 0
    for choice in itertools.product(range(len(per_pos)), repeat=Q):
        j = sum(1 for s in choice if per_pos[s][1] == -1)
        psi = {}
        for u in B:
            val = 1
            for s, idx in enumerate(choice):
                val *= per_pos[idx][0].get(u[s], 0)
                if val == 0:
                    break
            if val:
                psi[u] = val
        image = {}
        for w, c in psi.items():
            for u, cc in M0_cols[w].items():
                image[u] = image.get(u, 0) + c * cc
        image = {k: v for k, v in image.items() if v != 0}
        expected = (-1) ** j * _poly_eval(f_coefs, j)
        if image != {k: expected * v for k, v in psi.items()}:
            raise AssertionError("M0 does not act on the eigenbasis as V0 f(W0) should")
        if j in I:
            fixed += 1
        else:
            C_rows.append(psi)
    c_rows = len(C_rows)

    # Step 3: rows P_u(T) from (2.16); check the three properties.
    P_aux = {}
    for u in B:
        row = {}
        for w in B:
            m0 = M0_rows[u].get(w, 0)
            const = (1 if w == u else 0) - m0
            coefs = [const] if const else []
            if m0:
                E = [s + 1 for s in range(Q) if w[s][0] != u[s][1]]          # E_{u,w}
                if len(E) > R - 1:
                    raise AssertionError("support condition |E_{u,w}| <= R - 1 of (2.14) fails")
                phi = [Fraction(1)]
                for s in E:
                    phi = _pmul(phi, [Fraction(1), Fraction(-1, s)])          # (1 - T/s)
                coefs = _padd(coefs, [m0 * c for c in phi])
            while coefs and coefs[-1] == 0:
                coefs.pop()
            if coefs:
                row[w] = coefs
        if {w: cs[0] for w, cs in row.items() if cs[0] != 0} != {u: 1}:
            raise AssertionError("P_u(0) != e_u")
        if any(len(cs) - 1 > R - 1 for cs in row.values()):
            raise AssertionError("deg P_u > R - 1")
        P_aux[u] = row
        # (I - M0) row lies in the row space of C (first condition of (2.14))
        im = {w: (1 if w == u else 0) - M0_rows[u].get(w, 0) for w in B}
        im = {w: c for w, c in im.items() if c != 0}
        if _solve_in_row_space(im, C_rows) is None:
            raise AssertionError("row of I - M0 outside the row space of C")

    # Step 4: schedules a_v, d_e in [r]^Q, labels u(e, v), copies.
    scheds = list(itertools.product(range(r), repeat=Q))
    order = list(range(n))
    rng.shuffle(order)
    a = [None] * n
    for i, v in enumerate(order):
        a[v] = scheds[i % len(scheds)]
    d = []
    mult: dict = {}
    for e, c in enumerate(occ):
        best = None
        for cand in scheds:
            trial = dict(mult)
            for v, _ in c:
                lab = tuple((a[v][s], cand[s]) for s in range(Q))
                trial[lab] = trial.get(lab, 0) + 1
            key = (max(trial.values()), sum(x * x for x in trial.values()), rng.random())
            if best is None or key < best[0]:
                best = (key, cand, trial)
        d.append(best[1])
        mult = best[2]
    h = max(mult.values()) if mult else 0
    label = [[tuple((a[v][s], d[e][s]) for s in range(Q)) for v, _ in c] for e, c in enumerate(occ)]
    copy = [[0] * len(c) for c in occ]
    seen: dict = {}
    fill = [dict() for _ in range(h)]
    for e, c in enumerate(occ):
        for k, (v, _) in enumerate(c):
            lab = label[e][k]
            j = seen.get(lab, 0)
            seen[lab] = j + 1
            copy[e][k] = j
            fill[j][lab] = v
    # Replacement forms over the formula's variables, eq. (2.18).
    P = []
    for e, c in enumerate(occ):
        Pe = []
        for k, (v, _) in enumerate(c):
            form: dict = {}
            for w, coefs in P_aux[label[e][k]].items():
                if w in fill[copy[e][k]]:
                    var = fill[copy[e][k]][w]
                    form[var] = _padd(form.get(var, []), coefs)
            form = {var: cs for var, cs in form.items() if any(x != 0 for x in cs)}
            if {var: cs[0] for var, cs in form.items() if cs[0] != 0} != {v: 1}:
                raise AssertionError("P_{e,v}(0) != e_v after substitution")
            Pe.append(form)
        P.append(Pe)
    # Designs at every point: S_{t,g}, sigma_e(t), and the group-split measurement rows.
    points = {}
    for t in range(1, Q + 1):
        group_of = [a[v][t - 1] for v in range(n)]
        groups = [[v for v in range(n) if group_of[v] == g] for g in range(r)]
        sigma = [d[e][t - 1] for e in range(len(occ))]
        rows = [[] for _ in range(r)]
        seen_rows = [set() for _ in range(r)]
        for j in range(h):
            for psi in C_rows:
                ell: dict = {}
                for w, cc in psi.items():
                    if w in fill[j]:
                        var = fill[j][w]
                        ell[var] = ell.get(var, 0) + cc
                for g in range(r):
                    part = {var: cc for var, cc in ell.items() if group_of[var] == g and cc != 0}
                    key = tuple(sorted(part.items()))
                    if part and key not in seen_rows[g]:
                        seen_rows[g].add(key)
                        rows[g].append(part)
        points[t] = DesignPoint(groups, group_of, sigma, rows)
    info = {"D": Dn, "c": c_rows, "fixed": fixed, "I": I, "f": f_coefs, "h": h, "copies": h,
            "labels_used": len(mult), "a": a, "d": d, "C_rows": C_rows}
    return ToyDesign("aux", r, Q, R - 1, occ, P, points, M, info)


def brute_force_H(cnf: CNF, D: ToyDesign, t: int):
    """H(t) = sum_x prod_e Phi_e(<P_{e,v}(t), x>), by enumeration."""
    n = cnf.n
    forms = [[{i: _poly_eval(cs, t) for i, cs in D.P[e][k].items()} for k in range(len(D.occ[e]))] for e in range(len(D.occ))]
    total = 0
    for x in range(1 << n):
        bits = [(x >> i) & 1 for i in range(n)]
        prod = 1
        for e, c in enumerate(D.occ):
            inputs = [sum((coef * bits[i] for i, coef in forms[e][k].items()), 0) for k in range(len(c))]
            prod *= clause_poly([s for _, s in c], inputs, 1)
            if prod == 0:
                break
        total += prod
    return total


@dataclass
class ToyCounters:
    tuples: int = 0
    assignment_evals: int = 0
    clause_evals: int = 0
    restricted_sums: int = 0


def _prepare_point(D: ToyDesign, t: int):
    """Everything the measurement algorithm needs at t: numeric forms and, per group, the
    list of (assignment, measurement vector, own-group input values) for every assignment."""
    pt = D.at(t)
    own, lam = D.numeric_forms(t)
    pre = []
    for g, S in enumerate(pt.groups):
        clauses_g = [e for e in range(len(D.occ)) if pt.sigma[e] == g]
        entries = []
        for bits in range(1 << len(S)):
            u = {S[i]: (bits >> i) & 1 for i in range(len(S))}
            y = tuple(sum(c * u[v] for v, c in row.items()) for row in pt.rows[g])
            ownvals = {(e, k): sum((c * u[i] for i, c in own[e][k].items()), 0) for e in clauses_g for k in range(len(D.occ[e]))}
            entries.append((u, y, ownvals))
        pre.append((clauses_g, entries))
    return pt, own, lam, pre


def measurement_algorithm(D: ToyDesign, t: int, one, corr: Optional[dict] = None,
                          fixed_one: frozenset = frozenset(), counters: Optional[ToyCounters] = None,
                          cache: Optional[dict] = None):
    """Identity (2.6): Z = sum_y prod_g B_g(y). `one` is the ring unit (1 or a Poly); `corr`
    maps an occurrence (e, k) to a ring element added to its input (the <Delta P, 1_U>
    correction of (2.7)); `fixed_one` restricts the sum to assignments with x_i = 1 on E.
    The fibre loop `for ownvals in fibers[g][y[g]]` is the inner sum of (2.6)."""
    import itertools
    if cache is not None and t in cache:
        pt, own, lam, pre = cache[t]
    else:
        pt, own, lam, pre = _prepare_point(D, t)
        if cache is not None:
            cache[t] = (pt, own, lam, pre)
    fibers = []
    for g, (clauses_g, entries) in enumerate(pre):
        fib: dict = {}
        for u, y, ownvals in entries:
            if any(u[v] == 0 for v in pt.groups[g] if v in fixed_one):
                continue
            fib.setdefault(y, []).append(ownvals)
        fibers.append(fib)
    zero = one * 0
    Z = zero
    if counters:
        counters.restricted_sums += 1
    for y in itertools.product(*[list(f.keys()) for f in fibers]):
        if counters:
            counters.tuples += 1
        prod = one
        for g, (clauses_g, _entries) in enumerate(pre):
            foreign = {}
            for e in clauses_g:
                for k in range(len(D.occ[e])):
                    val = 0
                    for hh, coefs in lam[e][k].items():
                        val += sum((c * yy for c, yy in zip(coefs, y[hh])), 0)
                    foreign[(e, k)] = val
            Bg = zero
            for ownvals in fibers[g][y[g]]:
                if counters:
                    counters.assignment_evals += 1
                Ag = one
                for e in clauses_g:
                    inputs = []
                    for k in range(len(D.occ[e])):
                        val = ownvals[(e, k)] + foreign[(e, k)]
                        if corr and (e, k) in corr:
                            val = val + corr[(e, k)]
                        inputs.append(val)
                    Ag = Ag * clause_poly([s for _, s in D.occ[e]], inputs, one)
                    if counters:
                        counters.clause_evals += 1
                Bg = Bg + Ag
            prod = prod * Bg
        Z = Z + prod
    return Z


def taylor_data(D: ToyDesign, n: int, t: int, q: int, counters: ToyCounters, cache: dict) -> list:
    """Identity (2.9): the coefficients of z^0..z^q in H(t + z), from restricted sums Z_{t,E,U}."""
    import itertools
    J = q + 1
    one = Poly.const(1, J)
    delta: dict = {}
    for e in range(len(D.occ)):
        for k in range(len(D.occ[e])):
            dd = {}
            for i, cs in D.P[e][k].items():
                sh = _taylor_shift(cs, t)
                sh[0] = 0                                  # Delta P(z) = P(t + z) - P(t)
                pz = Poly(sh, J)
                if not pz.is_zero():
                    dd[i] = pz
            delta[(e, k)] = dd
    H = Poly([], J)
    for size in range(0, min(q, n) + 1):
        for E in itertools.combinations(range(n), size):
            for usize in range(size + 1):
                for U in itertools.combinations(E, usize):
                    corr = {}
                    for key, dd in delta.items():
                        s = None
                        for i in U:
                            if i in dd:
                                s = dd[i] if s is None else s + dd[i]
                        if s is not None:
                            corr[key] = s
                    Z = measurement_algorithm(D, t, one, corr, frozenset(E), counters, cache)
                    H = H + Z * ((-1) ** (size - usize))
    return H.coeffs()


def hermite_interpolate(data: dict, deg: int) -> Fraction:
    """Solve for the coefficients of H (degree <= deg) from Taylor data h_{t,j} = [z^j] H(t+z)
    and return H(0). Exact Gaussian elimination over Q."""
    rows = []
    for (t, j), hv in data.items():
        rows.append([Fraction(math.comb(i, j) * t ** (i - j)) if i >= j else Fraction(0) for i in range(deg + 1)] + [Fraction(hv)])
    ncol = deg + 1
    piv_row = 0
    pivots = []
    for col in range(ncol):
        p = next((i for i in range(piv_row, len(rows)) if rows[i][col] != 0), None)
        if p is None:
            continue
        rows[piv_row], rows[p] = rows[p], rows[piv_row]
        pv = rows[piv_row][col]
        rows[piv_row] = [x / pv for x in rows[piv_row]]
        for i in range(len(rows)):
            if i != piv_row and rows[i][col] != 0:
                f = rows[i][col]
                rows[i] = [a - f * b for a, b in zip(rows[i], rows[piv_row])]
        pivots.append(col)
        piv_row += 1
    for i in range(piv_row, len(rows)):
        if rows[i][-1] != 0:
            raise ValueError("inconsistent Taylor data: the degree bound or the identities failed")
    if len(pivots) < ncol:
        raise ValueError("underdetermined: Q (q+1) must exceed deg H")
    return rows[0][-1]


def _fmt_row(row: dict) -> str:
    parts = []
    for v, c in sorted(row.items()):
        sign = "-" if c < 0 else "+"
        mag = "" if abs(c) == 1 else f"{abs(c)}"
        parts.append(f"{sign} {mag}x{v + 1}")
    s = " ".join(parts)
    return s[2:] if s.startswith("+ ") else s


def toy_count(cnf: CNF, r: int, kind: str, deform: bool, q_target: int, seed: int,
              Qaux: int = 2, Raux: int = 2, budget: int = 2000, force: bool = False) -> None:
    rng = random.Random(seed)
    n = cnf.n
    print(f"formula: n = {n}, m = {cnf.m}, M = {sum(len(c) for c in cnf.clauses)} literal occurrences, bounded occurrence verified")
    t0 = time.perf_counter()
    true_count = count_solutions(cnf)
    t_bf = time.perf_counter() - t0
    print(f"brute force: #F = {true_count}  ({1 << n} assignments x {cnf.m} clauses, {t_bf*1e3:.0f} ms)")

    if kind == "direct":
        D = build_direct_design(cnf, r, deform, rng)
        if D.degP == 0:
            Q, q = 1, 0
        else:
            Q = (D.M * D.degP) // (max(0, q_target) + 1) + 1
            q = (D.M * D.degP) // Q
        D.Q = Q
    elif kind == "collapse":
        D = build_collapse_design(cnf, r)
        Q, q = 1, D.M * D.degP
    elif kind == "aux":
        D = build_aux_design(cnf, r, Qaux, Raux, rng)
        Q, q = D.Q, (D.M * D.degP) // D.Q
    else:
        raise ValueError(kind)
    J = q + 1
    degH = D.M * D.degP

    print(f"design: {kind}, r = {D.r} groups, degP = {D.degP}, deg H <= M degP = {degH}")
    if kind == "aux":
        inf = D.info
        print(f"        Section 2.7 construction: Q = {D.Q}, R = {Raux} (deg f <= {Raux - 1}), auxiliary coordinates D = r^(2Q) = {inf['D']},")
        print(f"        f interpolates (-1)^j on I = {inf['I']}: f(j) = {' + '.join(f'({c}) j^{i}' for i, c in enumerate(inf['f']))}")
        print(f"        eigenvectors fixed by M0: {inf['fixed']} of {inf['D']}; measurement rows of C: c = {inf['c']} (each a +-1 vector on {2 ** D.Q if D.r == 2 else 'several'} coordinates)")
        print(f"        schedules: {inf['labels_used']} labels used by {D.M} occurrences, max multiplicity h = {inf['h']} -> {inf['copies']} copies")
        for t in range(1, D.Q + 1):
            pt = D.at(t)
            print(f"        t = {t}: groups {[[v + 1 for v in S] for S in pt.groups]}, clause -> group {[g + 1 for g in pt.sigma]}")
            for g in range(len(pt.groups)):
                print(f"               C_{{{t},{g + 1}}}: {len(pt.rows[g])} rows: " + "; ".join(_fmt_row(row) for row in pt.rows[g]))
        ex_e, ex_k = next(((e, k) for e in range(len(D.occ)) for k in range(len(D.occ[e])) if len(D.P[e][k]) > 1), (0, 0))
        print(f"        example form, clause {ex_e + 1} occurrence of x{D.occ[ex_e][ex_k][0] + 1}: " +
              " + ".join(f"({' + '.join(f'{c} T^{i}' for i, c in enumerate(cs))}) x{v + 1}" for v, cs in sorted(D.P[ex_e][ex_k].items())))
    else:
        pt = D.at(1)
        print(f"        groups {[[v + 1 for v in S] for S in pt.groups]}, clause -> group {[g + 1 for g in pt.sigma]}, rows per group {[len(R_) for R_ in pt.rows]}"
              + (f", boundary variables read directly: {D.info['boundary']}" if kind == "direct" else ", one all-ones row per group"))
    print(f"        Q = {Q} points, q = floor(M degP / Q) = {q}, J = q + 1 = {J}, Q J = {Q * J} > {degH}")

    counters = ToyCounters()
    cache: dict = {}
    print("identity (2.6) at the design points:")
    t_26 = time.perf_counter()
    Z = None
    values = {}
    tuples_per_point = []
    for t in range(1, Q + 1):
        c0 = counters.tuples
        Z = measurement_algorithm(D, t, 1, None, frozenset(), counters, cache)
        values[t] = Z
        tuples_per_point.append(counters.tuples - c0)
        Ht = brute_force_H(cnf, D, t)
        print(f"  t = {t:2d}: measurement algorithm Z = {str(Z):>14}  brute-force H(t) = {str(Ht):>14}  {'OK' if Z == Ht else 'FAIL'}")
        if Z != Ht:
            raise SystemExit("identity (2.6) failed")
    pt1 = D.at(1)
    box = 1
    for g in range(len(pt1.groups)):
        Dg = max((sum(abs(c) for c in row.values()) for row in pt1.rows[g]), default=0)
        box *= (2 * Dg + 1) ** len(pt1.rows[g])
    tup = str(tuples_per_point) if Q <= 4 else f"{min(tuples_per_point)}..{max(tuples_per_point)}"
    print(f"  measurement tuples attained {tup} per point (paper's box bound prod (2D+1)^c_g = {box} at t = 1); "
          f"assignment evaluations {counters.assignment_evals}, clause evaluations {counters.clause_evals}, "
          f"{(time.perf_counter() - t_26)*1e3:.0f} ms")

    pairs_per_point = sum(math.comb(n, s) * (1 << s) for s in range(min(q, n) + 1))
    group_cost = sum(1 << len(S) for S in pt1.groups)
    print(f"the three costs of (2.11) at this scale: Taylor pairs (E, U) per point {pairs_per_point} x measurement tuples "
          f"{max(tuples_per_point)} x assignments within groups {group_cost}  (brute force: 2^n = {1 << n})")

    if D.degP == 0:
        print("identity (2.9): degP = 0, so H is constant and H(0) = H(1); no Taylor data needed")
        H0 = Z
    elif q == 0:
        print(f"identity (2.9): q = 0, so the only restricted sum per point is Z_(t,{{}},{{}}) = H(t), already computed above")
        H0 = hermite_interpolate({(t, 0): values[t] for t in range(1, Q + 1)}, degH)
        print(f"plain interpolation over {Q} values for a polynomial of degree <= {degH}: H(0) = {H0}")
    elif pairs_per_point * Q > budget and not force:
        print(f"identity (2.9): q = {q} >= {'n' if q >= n else 'n/2'} here, so the Taylor stage would need {pairs_per_point} restricted sums per point "
              f"({pairs_per_point * Q} in all, over the budget of {budget}; pass --force to run it).")
        print("  H(0) not recovered on this run. The pipeline is certified above at every design point; the Taylor stage is the")
        print("  same code that recovers H(0) exactly with --design direct, where q < n/2.")
        H0 = None
    else:
        note = " (q >= n: exhaustive at this scale)" if q >= n else ""
        print(f"identity (2.9): per point, {sum(math.comb(n, s) for s in range(min(q, n) + 1))} sets E with |E| <= {q} and {pairs_per_point} pairs (E, U){note}; series modulo z^{J}")
        counters = ToyCounters()
        data = {}
        t_29 = time.perf_counter()
        for t in range(1, Q + 1):
            coeffs = taylor_data(D, n, t, q, counters, cache)
            for j, hv in enumerate(coeffs):
                data[(t, j)] = hv
            if Q <= 4 or t in (1, Q):
                print(f"  t = {t:2d}: Taylor coefficients {[str(c) for c in coeffs[:min(6, len(coeffs))]]}{' ...' if len(coeffs) > 6 else ''}")
        print(f"  restricted sums {counters.restricted_sums}, tuples {counters.tuples}, assignment evaluations {counters.assignment_evals}, "
              f"clause evaluations {counters.clause_evals}, {(time.perf_counter() - t_29)*1e3:.0f} ms")
        H0 = hermite_interpolate(data, degH)
        print(f"Hermite interpolation over {len(data)} conditions for a polynomial of degree <= {degH}: H(0) = {H0}")
    if H0 is not None:
        ok = H0 == true_count
        print(f"result: H(0) = {H0} vs brute-force #F = {true_count}  {'OK' if ok else 'FAIL'}")
        if not ok:
            raise SystemExit("H(0) != #F")
    print("cost note: this toy certifies the identities and shows their cost; nothing here beats 2^n. The paper's saving needs")
    print("           R << Q with r^(2Q) << n (its Section 2.6), which exists only for enormous n.")


def main(argv: Optional[list[str]] = None) -> int:
    p = argparse.ArgumentParser(description="toy of OpenAI Theorem 1.1: identities (2.6), (2.9) and interpolation, checked by brute force")
    p.add_argument("--n", type=int, default=10)
    p.add_argument("--m", type=int, default=10, help="clauses (at most n for the bounded-occurrence form)")
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("--dimacs", help="use this CNF instead of a random one (must be bounded occurrence)")
    p.add_argument("--r", type=int, default=2, help="number of groups")
    p.add_argument("--design", choices=["direct", "collapse", "aux"], default="direct",
                   help="direct: unit rows on boundary variables; collapse: one all-ones row per group; aux: the paper's Section 2.7 construction")
    p.add_argument("--no-deform", action="store_true", help="direct design without the own-group deformation (degP = 0)")
    p.add_argument("--q", type=int, default=1, help="Taylor order per point for the direct design (Q is derived from it; 0 = plain interpolation)")
    p.add_argument("--Q", type=int, default=2, help="aux design: number of evaluation points (D = r^(2Q) auxiliary coordinates)")
    p.add_argument("--R", type=int, default=2, help="aux design: degree parameter (deg f <= R - 1; needs R <= Q + 1)")
    p.add_argument("--budget", type=int, default=2000, help="largest number of restricted sums the Taylor stage may run")
    p.add_argument("--force", action="store_true", help="run the Taylor stage even over the budget")
    args = p.parse_args(argv)
    if args.dimacs:
        with open(args.dimacs) as fh:
            cnf = CNF.from_dimacs(fh.read())
    else:
        cnf = random_bounded_occurrence_3cnf(args.n, args.m, random.Random(args.seed))
    toy_count(cnf, args.r, args.design, not args.no_deform, args.q, args.seed,
              Qaux=args.Q, Raux=args.R, budget=args.budget, force=args.force)
    return 0


if __name__ == "__main__":
    sys.exit(main())
