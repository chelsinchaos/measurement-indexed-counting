# Notes on family 378 — *The Exponential Time Hypothesis is False* (preprint of 25 October 2026)

Contributed as community commentary. These notes do not bear on the correctness of the theorem, whose formal verification is the Lean development `lean/docs/378.md`; they add an executable, inspectable model of the algorithm at toy scale, two structural remarks about the construction of Section 2.7, and two expository suggestions. Everything numeric below is reproducible from the linked repository, where it is checked by three independent runs.

## 1. An executable model of Sections 2.3–2.7

`checks/measurement_toy.py` (Python 3.10+, standard library only, exact arithmetic) implements, for a bounded-occurrence 3-CNF:

- the clause polynomials Φ_e and replacement linear forms P_{e,v}(T) (Section 2.3);
- measurement designs in the sense of Definition 2.1, as per-point objects (groups S_{t,g}, clause assignments σ_e(t), integer rows C_{t,g}); condition (4) is not assumed but verified at every point by solving (2.5) for the coefficients λ_{e,v,h};
- the measurement-tuple algorithm and identity (2.6);
- the restricted sums Z_{t,E,U} of (2.8) over Q[z]/(z^J), the inclusion–exclusion identity (2.9), and Hermite interpolation of H from the Taylor data at t = 1, …, Q (Section 2.5);
- three families of designs, including the paper's own Section 2.7 construction at small parameters (auxiliary coordinates B = ([r]×[r])^Q, swap matrices, M_0 = V_0 f(W_0) with f interpolating (−1)^j on an interval of R consecutive integers, measurement rows from the unfixed eigenvectors, rows P_u(T) of (2.16), schedules and copies of Step 4, group-split matrices C_{t,g}). The construction is checked as it is built (action of M_0 on the eigenbasis, support condition of (2.14), P_u(0) = e_u, degree bound, rows of I − M_0 in the row space of C).

Every stage is compared with brute force. Example runs (`python3 checks/measurement_toy.py --n 12 --m 12 --q 1`; `--design aux`):

| run | design | Q, q | (2.6) certified | H(0) = #F |
|---|---|---|---|---|
| n = 10, M = 27 | direct, deformed (deg P = 1) | 14, 1 | 14 of 14 points | 259 ✓ (Taylor, 294 restricted sums) |
| n = 12, M = 33 | direct, deformed | 17, 1 | 17 of 17 | 736 ✓ |
| n = 14, M = 39 | direct, deformed | 20, 1 | 20 of 20 | 3,108 ✓ (580 restricted sums, 40 Hermite conditions) |
| n = 6, M = 15 | Section 2.7 at r = 2, Q = 2, R = 2 | 2, 7 (exhaustive) | 2 of 2 | 34 ✓ |
| n = 12, M = 33 | Section 2.7 at r = 2, Q = 2, R = 2 | 2, 16 | 2 of 2 (H(1) = −12561/32, rational, equal to brute force) | not run (Taylor stage exhaustive here) |
| n = 12 | Section 2.7 at (r, Q, R) = (2, 3, 2) and (3, 2, 2) | 3 / 2 | all | not run |

At n = 12 with r = 2, Q = 2, R = 2 the construction yields D = 16 coordinates, f(j) = 1 − 2j, 15 of 16 eigenvectors fixed, one measurement row, 3 copies, and group-split rows that are ±1 differences of two variables (for instance x₃ − x₁₂ and x₇ − x₈ at t = 1); groups, assignments and rows all change between t = 1 and t = 2, and inputs at t ≠ 0 are rational (e.g. −½·T·x₃). A reader who wants to *see* the mechanism of Section 2.7 can print it from this run.

Independently of the toy, a second script expands H(T) symbolically by brute force over all 2^n assignments and confirms every recorded value H(t), every recorded Taylor coefficient and H(0), and re-tests the design conditions with a fresh rank routine.

## 2. Two structural remarks on the construction of Section 2.7

**(a) Within-group components.** Under the random schedules a_v, d_e ∈ [r]^Q of Step 4, consider at a point t the hypergraph in which a clause assigned to group g couples those of its variables that lie in S_{t,g} (the only coupling inside a group sum of (2.6); every other dependence enters through the fixed measurements). Exploring it from a variable, each of the variable's three clauses is assigned to g with probability 1/r, each of the other two variables of such a clause lies in g with probability 1/r, and a variable reached through a clause has two further clauses: the offspring mean is 4/r². The exploration is critical at r = 2 and subcritical for r ≥ 3. Measured on random bounded-occurrence 3-CNFs with m = n under uniform schedules (three seeds): the largest within-group component is 83–109 at n = 1,000, 196–229 at n = 4,000 and ≈ n^{2/3} at r = 2, and stays at ≤ 15 variables for r ≥ 3 up to n = 16,000. In the paper's regime r → ∞, so its group sums factorize over components of bounded size: the measurement y_g = C_{t,g} u is additive across components, and the group histogram over the measurement box is the convolution of the component histograms.

**(b) This does not change the exponent.** With the paper's choices (2.19), the three terms of its exponent n[(R/Q) log(2Q/R) + τ log(2D+1) + 1/r] are, to leading order, 256 ν²/√λ (Taylor data), ν²/(128√λ) (measurements) and ν²/√λ (group enumeration), with λ = ⌈log₂ n⌉, ν = ⌈log₂ λ⌉. Replacing the group enumeration of Step 2 by a component product would therefore change the exponent by about one part in 256. The remark is structural, not algorithmic.

## 3. Two expository suggestions

1. Definition 2.1 is stated at a fixed t ("Throughout this subsection, write S_g = S_{t,g} and C_g = C_{t,g}"). It may help readers to say explicitly, at the definition, that groups, clause assignments and measurement matrices are all allowed to depend on t, since a t-independent C_h cannot work: if the rows of C_h did not depend on t and Q > R, then for an occurrence with v ∈ S_h the vectors P_{e,v}(t)|_{S_h} = e_v + Σ_j t^j w_j|_{S_h}, t = 1, …, Q, would span e_v itself by interpolation, forcing C_h to read every boundary variable directly. (Our toy's "direct" design is exactly that degenerate case, and it compresses nothing.)
2. The degenerate case R = Q + 1 of Section 2.7 (f fixes every eigenvector, c = 0, no measurements) is a clean illustration of the sentence in Section 2.7 that interpolating a different form at each point "would generally cost degree Q"; it runs in the toy as `--design aux --Q 3 --R 4`.

## 4. Where the rest is

Full report (methods, tables, the dense-vs-sparse recombination experiment, the semiring-cache experiment), code, dataset and verification logs: https://github.com/chelsinchaos/measurement-indexed-counting. Author: Chelsea A. McElveen. Implementation and drafting assisted by Claude (Anthropic); the author is responsible for the content.
