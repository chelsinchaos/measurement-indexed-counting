# PR: family 378 — executable toy model, structural notes, expository suggestions

**Folder:** `preprints/The-Exponential-Time-Hypothesis-is-False-October-25-2026/`

**Adds**
- `NOTES.md` — commentary: an executable model of Sections 2.3–2.7 checked against brute force at n ≤ 14; two structural remarks on the Section 2.7 construction (within-group components are critical at r = 2 and subcritical for r ≥ 3; the exponent accounting showing this does not change the result); two expository suggestions (making the t-dependence of Definition 2.1 explicit, with the one-line interpolation argument for why it is necessary; the degenerate case R = Q + 1).
- `checks/measurement_toy.py` — the model (Python 3.10+, standard library only, exact arithmetic); `python3 checks/measurement_toy.py --n 12 --m 12 --design aux` prints the Section 2.7 design at r = 2, Q = 2, R = 2 and certifies (2.6) at both points.
- `checks/README.md` — how to run, what is checked, pointer to the full repository with dataset and three verification runs.

**Does not change** the TeX source or the manuscript's claims. No error in the manuscript is asserted.

**Why it belongs here.** The README invites exposition improvements and commentary. A toy in which the measurement rows (e.g. x₃ − x₁₂), the rational inputs (−½·T·x₃) and the intermediate values can be printed is a comprehension aid that complements the Lean formalization.

**Review notes.** All numbers in NOTES.md are reproduced by the script; the independent re-computation and the dataset with SHA-256 manifest are in the linked repository. Contributed by Chelsea A. McElveen; implementation and drafting assisted by Claude (Anthropic).
