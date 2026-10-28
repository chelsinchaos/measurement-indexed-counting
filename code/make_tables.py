#!/usr/bin/env python3
"""
make_tables.py - produce every table (LaTeX) and figure (PDF) of the paper from the dataset.
Nothing numeric in the paper is typed by hand: the paper \\input{}s the files written here.
"""
from __future__ import annotations

import json
import math
import os
import sys
from fractions import Fraction

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
DS = os.path.join(HERE, "..", "dataset")
OUT = os.path.join(HERE, "..", "paper", "tables")
FIG = os.path.join(HERE, "..", "paper", "figures")
os.makedirs(OUT, exist_ok=True)
os.makedirs(FIG, exist_ok=True)


def load(name):
    return json.load(open(os.path.join(DS, "results", name)))["data"]


def w(name, text):
    with open(os.path.join(OUT, name), "w") as fh:
        fh.write(text)


def num(x):
    return f"{x:,}"


def ms(t):
    return f"{t*1e3:,.0f}" if t * 1e3 >= 100 else f"{t*1e3:.1f}"


def secs(t):
    return f"{t:,.1f}"


def tex_escape(s):
    return s.replace("_", r"\_").replace("%", r"\%").replace("&", r"\&")


macros = []


def macro(name, value):
    macros.append(f"\\newcommand{{\\{name}}}{{{value}}}")


# ---------------------------------------------------------------------------
# E1: the executable model
# ---------------------------------------------------------------------------
E1 = load("E1_toy.json")
runs = {r["name"]: r for r in E1["runs"]}
order = ["toy_direct_n8_q1", "toy_direct_n10_q1", "toy_direct_n12_q1", "toy_direct_n14_q1", "toy_direct_n14_q0", "toy_direct_n12_nodeform",
         "toy_aux_n6_Q2R2", "toy_aux_n12_Q2R2", "toy_aux_n12_Q3R2", "toy_aux_n12_r3_Q2R2", "toy_aux_n8_Q3R4", "toy_collapse_n7"]
label = {
    "toy_direct_n8_q1": "direct, deformed", "toy_direct_n10_q1": "direct, deformed", "toy_direct_n12_q1": "direct, deformed", "toy_direct_n14_q1": "direct, deformed",
    "toy_direct_n14_q0": "direct, deformed, values only", "toy_direct_n12_nodeform": "direct, undeformed ($\\deg P = 0$)",
    "toy_aux_n6_Q2R2": "aux (\\S2.7), $r{=}2,Q{=}2,R{=}2$", "toy_aux_n12_Q2R2": "aux (\\S2.7), $r{=}2,Q{=}2,R{=}2$",
    "toy_aux_n12_Q3R2": "aux (\\S2.7), $r{=}2,Q{=}3,R{=}2$", "toy_aux_n12_r3_Q2R2": "aux (\\S2.7), $r{=}3,Q{=}2,R{=}2$",
    "toy_aux_n8_Q3R4": "aux (\\S2.7), $r{=}2,Q{=}3,R{=}4{=}Q{+}1$", "toy_collapse_n7": "collapse",
}
rows = []
for name in order:
    if name not in runs:
        continue
    r = runs[name]
    ts = r.get("taylor_stage", {})
    stage = ts.get("stage", "")
    if stage == "taylor_2_9":
        rec = f"{r['H0']} (Taylor, {num(ts['restricted_sums'])} restricted sums)"
    elif stage == "plain_interpolation":
        rec = f"{r['H0']} (values only)"
    elif stage == "constant":
        rec = f"{r['H0']} (constant $H$)"
    else:
        rec = f"not run ({num(ts.get('restricted_sums_needed', 0))} sums)"
    okmark = "\\checkmark" if r.get("H0_equals_count", False) else "--"
    rows.append(f"{r['n']} & {r['m']} & {r['M']} & {label[name]} & {r['Q']} & {r['q']} & {len(r['identity_2_6']['points'])} of {r['Q']} & "
                f"{r['brute_force_count']} & {rec} & {okmark} & {secs(r['t_total_s'])} \\\\")
w("E1_runs.tex", "\n".join(rows) + "\n")

# costs table
rows = []
for name in order:
    if name not in runs:
        continue
    r = runs[name]
    c = r["costs_2_11"]
    rows.append(f"{r['n']} & {label[name]} & {num(c['taylor_pairs_per_point'])} & {num(c['tuples_max'])} & {num(c['box_bound_t1'])} & "
                f"{num(c['group_assignments'])} & {num(c['two_to_n'])} \\\\")
w("E1_costs.tex", "\n".join(rows) + "\n")

# the n = 12 aux example design (rows and example form) for the text
ex = json.load(open(os.path.join(DS, "designs", "toy_aux_n12_Q2R2.json")))


def fmt_row(row):
    terms = []
    for v, c in sorted(row.items(), key=lambda kv: int(kv[0])):
        sign = "-" if c < 0 else "+"
        mag = "" if abs(c) == 1 else str(abs(c))
        terms.append(f"{sign} {mag}x_{{{v}}}")
    s_ = " ".join(terms)
    return s_[2:] if s_.startswith("+ ") else s_.replace("- ", "-", 1) if s_.startswith("- ") else s_


lines = ["\\begin{itemize}[leftmargin=1.4em,itemsep=1pt,topsep=2pt]"]
for t, pt in ex["points"].items():
    groups = ", ".join("\\{" + ",".join(str(v) for v in S) + "\\}" for S in pt["groups"])
    sigma = ",".join(str(g) for g in pt["sigma"])
    rows_txt = "; ".join(f"$C_{{{t},{g+1}}}$: " + ", ".join("$" + fmt_row(row) + "$" for row in rows_g) for g, rows_g in enumerate(pt["rows"]))
    lines.append(f"\\item $t = {t}$: groups ${groups}$; clause $\\to$ group $({sigma})$; measurement rows {rows_txt}.")
lines.append("\\end{itemize}")
w("E1_aux_example.tex", "\n".join(lines) + "\n")


def fmt_poly(cs):
    """coefficient strings (low to high in T) -> LaTeX polynomial with proper signs and fractions"""
    from fractions import Fraction
    parts = []
    for i, c in enumerate(cs):
        f = Fraction(c)
        if f == 0:
            continue
        sign = "-" if f < 0 else "+"
        a = abs(f)
        if a.denominator == 1:
            mag = "" if (a == 1 and i > 0) else str(a.numerator)
        else:
            mag = f"\\tfrac{{{a.numerator}}}{{{a.denominator}}}"
        term = mag + ("" if i == 0 else ("T" if i == 1 else f"T^{{{i}}}"))
        parts.append((sign, term))
    out = ""
    for j, (sign, term) in enumerate(parts):
        out += (("-" if sign == "-" else "") if j == 0 else f" {sign} ") + term
    return out or "0"


macro("auxTwelveD", str(ex["info"]["D"]))
macro("auxTwelveC", str(ex["info"]["c"]))
macro("auxTwelveH", str(ex["info"]["h"]))
macro("auxTwelveFixed", str(ex["info"]["fixed"]))
macro("auxTwelveLabels", str(ex["info"]["labels_used"]))
# an example form with a rational coefficient: the one touching the most variables
best = None
for e, fe in enumerate(ex["forms"]):
    for k, form in enumerate(fe):
        if any(isinstance(c, str) and "/" in c for cs in form.values() for c in cs):
            if best is None or len(form) > len(best[2]):
                best = (e, k, form)
e, k, form = best
parts = []
for v, cs in sorted(form.items(), key=lambda kv: int(kv[0])):
    poly = fmt_poly(cs)
    if poly == "1":
        parts.append(f"x_{{{v}}}")
    elif poly == "-1":
        parts.append(f"-x_{{{v}}}")
    else:
        parts.append(f"({poly})\\,x_{{{v}}}" if (" + " in poly or " - " in poly) else f"{poly}\\,x_{{{v}}}")
macro("auxTwelveExampleForm", " + ".join(parts).replace("+ -", "- "))
macro("auxTwelveExampleClause", str(e + 1))
own_var = [v for v, cs in form.items() if Fraction(cs[0]) == 1][0] if any(Fraction(cs[0]) == 1 for cs in form.values()) else "?"
macro("auxTwelveExampleVar", str(own_var))
p1 = runs["toy_aux_n12_Q2R2"]["identity_2_6"]["points"][0]
macro("auxTwelveHone", p1["Z"])
macro("auxTwelveHtwo", runs["toy_aux_n12_Q2R2"]["identity_2_6"]["points"][1]["Z"])
macro("auxTwelveCount", str(runs["toy_aux_n12_Q2R2"]["brute_force_count"]))
r14 = runs.get("toy_direct_n14_q1")
if r14:
    macro("directFourteenSums", num(r14["taylor_stage"]["restricted_sums"]))
    macro("directFourteenClauseEvals", num(r14["taylor_stage"]["clause_evals"]))
    macro("directFourteenMinutes", f"{r14['t_total_s']/60:.1f}")
    macro("directFourteenCount", str(r14["brute_force_count"]))
macro("toyRunsCount", str(len(E1["runs"])))
macro("auxSixQ", str(runs["toy_aux_n6_Q2R2"]["q"]))
macro("auxSixM", str(runs["toy_aux_n6_Q2R2"]["M"]))
macro("auxSixSums", num(runs["toy_aux_n6_Q2R2"]["taylor_stage"]["restricted_sums"]))
macro("toyPointsCertified", num(sum(len(r["identity_2_6"]["points"]) for r in E1["runs"])))

# ---------------------------------------------------------------------------
# E2: components
# ---------------------------------------------------------------------------
E2 = load("E2_components.json")
ns = sorted({r["n"] for r in E2})
rs = (2, 3, 4, 5, 6, 8, 12)
rows = []
for r in rs:
    cells = []
    for n in ns:
        vals = [x["largest"] for x in E2 if x["n"] == n and x["r"] == r]
        cells.append(f"{min(vals)}--{max(vals)}" if min(vals) != max(vals) else f"{min(vals)}")
    mu = 4 / r / r
    rows.append(f"{r} & {3/r:.2f} & {mu:.2f} & " + " & ".join(cells) + " \\\\")
w("E2_components.tex", "\n".join(rows) + "\n")
w("E2_header.tex", " & ".join(f"$n = {num(n)}$ (group ${num(n)}/r$)" for n in ns) + "\n")
# n^{2/3} comparison at r = 2
cmp = []
for n in ns:
    vals = [x["largest"] for x in E2 if x["n"] == n and x["r"] == 2]
    cmp.append(f"$n = {num(n)}$: largest {min(vals)}--{max(vals)}, $n^{{2/3}} = {n ** (2/3):.0f}$")
macro("compTwoThirds", "; ".join(cmp))

# figure: largest component vs n (r = 2, 3, 4) with n^{2/3}
fig, ax = plt.subplots(figsize=(4.6, 3.2))
for r, mk in ((2, "o"), (3, "s"), (4, "^"), (6, "d")):
    xs, ys = [], []
    for n in ns:
        for x in E2:
            if x["n"] == n and x["r"] == r:
                xs.append(n); ys.append(x["largest"])
    ax.scatter(xs, ys, marker=mk, label=f"$r = {r}$", s=22)
nn = [min(ns) / 1.2, max(ns) * 1.2]
ax.plot(nn, [v ** (2 / 3) for v in nn], "k--", lw=0.8, label="$n^{2/3}$")
ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xlabel("$n$ (variables; $m = n$ clauses)"); ax.set_ylabel("largest within-group component")
ax.legend(fontsize=7, frameon=False)
fig.tight_layout(); fig.savefig(os.path.join(FIG, "components.pdf")); plt.close(fig)

# ---------------------------------------------------------------------------
# E3: fibre loop
# ---------------------------------------------------------------------------
E3 = load("E3_fibre.json")


def e3_rows(recs):
    rows = []
    for r in recs:
        m = r["methods"]
        tb = ms(m["bucket"]["t_s"]) if "bucket" in m else "--"
        td = ms(m["dict"]["t_s"]) if "dict" in m else "--"
        tp = ms(m["packed"]["t_s"])
        rb = f"{m['bucket']['t_s']/m['packed']['t_s']:.2f}" if "bucket" in m else "--"
        rd = f"{m['dict']['t_s']/m['packed']['t_s']:.2f}" if "dict" in m else "--"
        comp = r["components"]
        dens = 100 * r["attained"] / r["lanes"]
        rows.append(f"{r['k']} & {r['rows']} & {num(r['lanes'])} & {comp['sizes'][0]} & {num(comp['enumerated'])} & {num(r['attained'])} ({dens:.0f}\\%) & "
                    f"{tb} & {td} & {tp} & {rb} & {rd} \\\\")
    return rows


w("E3_design.tex", "\n".join(e3_rows(E3["design_r2"])) + "\n")
syn = E3["synthetic"]
for r, rows_ in ((4, 8), (3, 8), (4, 10)):
    recs = [x for x in syn if x["r"] == r and x["synthetic_rows"] == rows_]
    w(f"E3_synth_r{r}_rows{rows_}.tex", "\n".join(e3_rows(recs)) + "\n")
allrecs = E3["design_r2"] + syn
ratios_d = [x["methods"]["dict"]["t_s"] / x["methods"]["packed"]["t_s"] for x in allrecs if "dict" in x["methods"]]
ratios_b = [x["methods"]["bucket"]["t_s"] / x["methods"]["packed"]["t_s"] for x in allrecs if "bucket" in x["methods"]]
macro("fibreInstances", str(len(allrecs)))
macro("fibreDictRatioMin", f"{min(ratios_d):.2f}")
macro("fibreDictRatioMax", f"{max(ratios_d):.2f}")
macro("fibreBucketRatioMin", f"{min(ratios_b):.2f}")
macro("fibreBucketRatioMax", f"{max(ratios_b):.1f}")
dens = [100 * x["attained"] / x["lanes"] for x in syn]
macro("fibreDensityMin", f"{min(dens):.0f}")
macro("fibreDensityMax", f"{max(dens):.0f}")
big = [x for x in syn if x["r"] >= 3 and x["k"] >= 20 and "bucket" in x["methods"]]
macro("fibreSubcritMin", f"{min(x['methods']['bucket']['t_s']/x['methods']['packed']['t_s'] for x in big):.1f}")
macro("fibreSubcritMax", f"{max(x['methods']['bucket']['t_s']/x['methods']['packed']['t_s'] for x in big):.1f}")
exact_checks = sum(1 for x in allrecs if x.get("exact_rational_check"))
macro("fibrePackedOverDictMax", f"{max(x['methods']['packed']['t_s'] / x['methods']['dict']['t_s'] for x in allrecs if 'dict' in x['methods']):.0f}")
def best_component(x):
    m = x["methods"]
    return min(m[k]["t_s"] for k in ("dict", "packed") if k in m)
sub = [x["methods"]["bucket"]["t_s"] / best_component(x) for x in syn if x["r"] >= 3 and x["k"] >= 20 and "bucket" in x["methods"]]
macro("fibreBestSubcritMin", f"{min(sub):.1f}")
macro("fibreBestSubcritMax", f"{max(sub):.1f}")
rtwo = [x["methods"]["bucket"]["t_s"] / best_component(x) for x in E3["design_r2"] if "bucket" in x["methods"]]
macro("fibreBestRtwoMax", f"{max(rtwo):.1f}")
macro("fibrePackedFasterCount", str(sum(1 for x in allrecs if "dict" in x["methods"] and x["methods"]["dict"]["t_s"] > x["methods"]["packed"]["t_s"])))
macro("fibreExactChecks", str(exact_checks))

# figure: time vs k, two panels
fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.0), sharey=True)
for ax, recs, title in ((axes[0], E3["design_r2"], "real Section 2.7 design, $r = 2$ (one giant component)"),
                        (axes[1], [x for x in syn if x["r"] == 4 and x["synthetic_rows"] == 8], "synthetic rows, $r = 4$ (components of at most two)")):
    for meth, mk in (("bucket", "o"), ("dict", "s"), ("packed", "^")):
        xs = [x["k"] for x in recs if meth in x["methods"]]
        ys = [x["methods"][meth]["t_s"] for x in recs if meth in x["methods"]]
        ax.plot(xs, ys, marker=mk, ms=4, lw=1, label=meth)
    ax.set_yscale("log"); ax.set_xlabel("group size $k$"); ax.set_title(title, fontsize=8)
    ax.grid(alpha=0.3)
axes[0].set_ylabel("seconds (one group, one foreign tuple)")
axes[0].legend(fontsize=7, frameon=False)
fig.tight_layout(); fig.savefig(os.path.join(FIG, "fibre.pdf")); plt.close(fig)

# ---------------------------------------------------------------------------
# E4: constrained counting
# ---------------------------------------------------------------------------
E4 = load("E4_constrained.json")
rows = []
for g in E4["grids"]:
    W, H = g["grid"]
    encs = "; ".join(f"$b{{=}}{e['b']}$: {e['aux']} aux, {ms(e['t_s'])} ms, {num(e['nodes'])} nodes, comp.\\ {e['max_component']} ({e['ratio_to_histogram']:.1f}$\\times$)" for e in g["encodings"]) or "--"
    oe = "= OEIS" if g["total_equals_oeis"] else ("" if g["total_equals_oeis"] is None else "DIFFERS")
    rows.append(f"${W}\\times{H}$ & {g['n']} & {len(g['histogram'])} & {ms(g['t_s'])} & {num(g['nodes'])} & {num(g['cache_entries'])}/{num(g['cache_hits'])} & {num(g['total'])} {oe} & {encs} \\\\")
w("E4_grids.tex", "\n".join(rows) + "\n")
rows = []
for r in E4["random"]:
    hs, ps = r["hist_stats"], r["plain_stats"]
    encs = "; ".join(f"$b{{=}}{e['b']}$: {e['aux']} aux, {ms(e['t_s'])} ms, {num(e['nodes'])} nodes, comp.\\ {e['max_component']} ({e['ratio_to_histogram']:.1f}$\\times$)" for e in r["encodings"])
    wt = "weights $w_1 \\in \\{1,2,3\\}$" if r["weighted"] else "unweighted"
    rows.append(f"{r['n']} & {r['m']} & {wt} & {len(r['histogram'])} & {ms(r['t_hist_s'])} & {num(hs['nodes'])} & {ms(r['t_plain_s'])} & {r['t_hist_s']/r['t_plain_s']:.2f} & {encs} \\\\")
w("E4_random.tex", "\n".join(rows) + "\n")
mr = {x["kind"]: x for x in E4["multirow"]}
g66 = next(g for g in E4["grids"] if g["grid"] == [6, 6])
macro("multirowPB", ms(mr["pseudo_boolean"]["t_s"]))
macro("multirowTwo", ms(mr["two_rows"]["t_s"]))
macro("multirowPBentries", str(len(mr["pseudo_boolean"]["histogram"])))
macro("multirowTwoEntries", str(len(mr["two_rows"]["histogram"])))
macro("gridSixSixMs", ms(g66["t_s"]))
all_enc = [e for g in E4["grids"] for e in g["encodings"]] + [e for r in E4["random"] for e in r["encodings"]]
macro("encRatioMin", f"{min(e['ratio_to_histogram'] for e in all_enc):.1f}")
macro("encRatioMax", num(round(max(e["ratio_to_histogram"] for e in all_enc))))
macro("encRatioGridMax", f"{max(e['ratio_to_histogram'] for g in E4['grids'] for e in g['encodings']):.0f}")
macro("encRatioRandomMax", num(round(max(e["ratio_to_histogram"] for r in E4["random"] for e in r["encodings"]))))
hr = [r["t_hist_s"] / r["t_plain_s"] for r in E4["random"]]
macro("histOverPlainMin", f"{min(hr):.1f}")
macro("histOverPlainMax", f"{max(hr):.1f}")
r40 = next((r for r in E4["random"] if r["n"] == 40), None)
if r40:
    e20 = next((e for e in r40["encodings"] if e["b"] == 20), None)
    if e20:
        macro("encFortyTwentySeconds", f"{e20['t_s']:.0f}")
        macro("encFortyTwentyNodes", num(e20["nodes"]))
        macro("encFortyTwentyComp", str(e20["max_component"]))
        macro("encFortyTwentyRatio", num(round(e20["ratio_to_histogram"])))
    macro("histFortyMs", ms(r40["t_hist_s"]))
    macro("histFortyEntries", str(len(r40["histogram"])))
g88 = next(g for g in E4["grids"] if g["grid"] == [8, 8])
macro("gridEightTotal", num(g88["total"]))
macro("gridEightMs", ms(g88["t_s"]))
macro("gridEightSizes", str(len(g88["histogram"])))
g66e = next(e for e in g66["encodings"] if e["b"] == g66["mode"])
macro("gridSixEncComp", str(g66e["max_component"]))
macro("gridSixEncRatio", f"{g66e['ratio_to_histogram']:.0f}")

# figure: encoding ratio vs b
fig, ax = plt.subplots(figsize=(4.6, 3.0))
for g in E4["grids"]:
    if g["encodings"]:
        ax.plot([e["b"] for e in g["encodings"]], [e["ratio_to_histogram"] for e in g["encodings"]], marker="o", ms=4, lw=1, label=f"grid ${g['grid'][0]}\\times{g['grid'][1]}$")
for r in E4["random"]:
    ax.plot([e["b"] for e in r["encodings"]], [e["ratio_to_histogram"] for e in r["encodings"]], marker="s", ms=4, lw=1, ls="--", label=f"random $n={r['n']}$")
ax.set_yscale("log"); ax.set_xlabel("constraint value $b$"); ax.set_ylabel("time of one encoded $b$ / whole histogram")
ax.axhline(1, color="k", lw=0.6)
ax.legend(fontsize=7, frameon=False); ax.grid(alpha=0.3)
fig.tight_layout(); fig.savefig(os.path.join(FIG, "encoding.pdf")); plt.close(fig)

# ---------------------------------------------------------------------------
# exponent terms under the paper's parameter choices (2.19), at finite lambda
# ---------------------------------------------------------------------------
def exponent_terms(lam):
    nu = math.ceil(math.log2(lam))
    r = math.floor(math.sqrt(lam) / nu ** 2)
    Q = math.floor(lam / (128 * nu))
    R = 4 * math.ceil(math.sqrt(lam))
    if r < 1 or Q < 1:
        return None
    T1 = (R / Q) * math.log2(2 * Q / R)
    T2 = (1 / (r * lam)) * (2 * Q * math.log2(r) + 1)       # tau * log2(2D+1), D = r^{2Q}
    T3 = 1 / r
    return T1, T2, T3
first_pos = next(k for k in range(10, 80) if exponent_terms(2 ** k) and exponent_terms(2 ** k)[0] > 0)
macro("expFirstPositive", str(first_pos))
# exact threshold: smallest integer lambda with r, Q >= 1 and T1 > 0 (T1 > 0 iff 2Q > R), by bisection
def t1_positive(lam):
    t = exponent_terms(lam)
    return t is not None and t[0] > 0
lo, hi = 2 ** (first_pos - 1), 2 ** first_pos
while hi - lo > 1:
    mid = (lo + hi) // 2
    if t1_positive(mid):
        hi = mid
    else:
        lo = mid
macro("expThresholdLambda", str(hi).replace(",", "")[:2] + "{,}" + str(hi)[2:5] + "{,}" + str(hi)[5:])
macro("expThresholdLog", f"{math.log2(hi):.1f}")
for k, tag in ((26, "A"), (32, "B"), (64, "C"), (128, "D"), (512, "E")):      # LaTeX macro names cannot carry digits
    T1, T2, T3 = exponent_terms(2 ** k)
    macro(f"expShare{tag}", f"{100 * T3 / (T1 + T2 + T3):.2g}" if k >= 64 else f"{100 * T3 / (T1 + T2 + T3):.1f}")
    macro(f"expRatio{tag}", f"{T1 / T3:.0f}")
macro("expShareLimit", f"{100 / (256 + 1 + 1 / 128):.2f}")

# ---------------------------------------------------------------------------
# independent checks + manifest
# ---------------------------------------------------------------------------
ic_path = os.path.join(DS, "results", "independent_check.json")
if os.path.exists(ic_path):
    ic = json.load(open(ic_path))
    macro("indepChecks", str(ic["checks"]))
    macro("indepPassed", str(ic["passed"]))
man = json.load(open(os.path.join(DS, "MANIFEST.json")))
macro("datasetFiles", str(len(man["files"])))
env = man["environment"]
macro("envPython", env["python"].split()[0])
macro("envNumpy", env["numpy"])
macro("envGmpy", env["gmpy2"])
macro("envCpus", str(env["cpu_count"]))
macro("envPlatform", tex_escape(env["platform"]))
w("macros.tex", "\n".join(macros) + "\n")
print(f"wrote {len(os.listdir(OUT))} table files and {len(os.listdir(FIG))} figures; {len(macros)} macros")
