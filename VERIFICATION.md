# VERIFICATION.md — independent audit (verification run 3)

Auditor: an agent with no access to the development history, working from the repository as found on 2026-10-28 (dataset generated 14:41–15:13 UTC; audit 15:17–15:30 UTC). No repository file was edited; the one file `verify_all.py` rewrites was restored byte-for-byte.

## Summary

**PASS.** Every number the paper takes from the dataset is what the dataset says, and the dataset is what the code produces: quick regeneration gave 0 differences in deterministic quantities, the independent re-implementation passed 74 of 74 checks, the manifest matched on all 64 files, all 12 table/macro files regenerate byte-identically, and 45 macros recomputed with my own code agree with `macros.tex`. Grid totals equal OEIS A006506 (also reproduced by my own transfer DP). The problems are in prose, not data: Table 1 omits one of the 12 runs the text says it lists; three ratio ranges are worded against their own values (a "win" of 0.4×, "0.6× more", "packed over dict" for a dict/packed ratio); one sentence says q = M where the table says q = 7, M = 15; and the Appendix B ratio 256 : 1 : 1/128 is algebraically right but unreachable at any writable n (the paper's parameters give r = 0 at n = 10^100 and 10^1000).

## Run 1/2 re-execution results

`cd code && python3 verify_all.py` (log in scratchpad `audit_verify_all.txt`); wall time **6 m 14 s**, not the advertised 15–20 min.

| Stage | Result |
|---|---|
| Manifest | `manifest: 64 files, 0 problems` |
| Run 1, `run_all.py --quick` | `regeneration: 0 differences in deterministic quantities` |
| Run 2, `independent_check.py` | `74 of 74 independent checks passed` (E1 12 runs; E2 42 records; E3 7 histograms k ≤ 16; E4 5 grids, 15 encodings, 2 multirow, 3 random) |
| Final line | `VERIFY_ALL: PASS` |

Beyond the script: all 38 regenerated `.cnf` and 11 `designs/*.json` are byte-identical to the released ones (`verify_all.py` compares results JSON only). Quick mode does **not** regenerate `toy_direct_n14_q1` (580 sums, 9.2 min), E2 at n = 16,000, the two k = 24 fibre instances (source of 103.6×), or the n = 40, b = 20 encoding (source of 1,083×); these rest on Run 2 alone (E1 symbolic check, E4 DPLL) or, for E2 n = 16,000 and the k = 24 timings, on no re-execution.

Side effect: `independent_check.py` overwrites `dataset/results/independent_check.json` in the released dataset (new `finished_utc` and timing strings; verdicts identical), so its SHA-256 stopped matching `MANIFEST.json`; `verify_all.py` silently exempts that file. I restored the original; all 64 hashes match again.

## Hand-typed numbers checked

| Location | Claim | Dataset | Verdict |
|---|---|---|---|
| Abstract | ≤ 15 of 1,333 variables at n = 4,000, r ≥ 3 | max 15 (r = 3); group sizes 1312–1369 | OK |
| Abstract, §1, §4.3 | n ≤ 14; design parameter list; "low double digits" at 16,000 | run list; max 19 | OK |
| §4.2, §4.4 | M = 33; 12D² = 3,072 at D = 16 | M = 33 | OK |
| §4.3 | Q = 20 points, 40 Hermite conditions, 9.2 min | Q = 20, J = 2, conditions 40, 552.8 s | OK |
| §4.3 | "n = 6, where Q = 2 forces q = M ≥ n" | q = 7, M = 15 (q = ⌊M/Q⌋ ≥ n makes the stage exhaustive: 2·3⁶ = 1,458 sums) | **wrong** |
| §5.1 | three seeds, m = n, 7 values of r | 63 = 3·7·3 records | OK |
| §6.1 | p = 2³¹ − 1; 9–10 rows per copy for r ≥ 3, Q = 2 | `PRIME`; c = 9 (r = 3, Q = 2), c = 10 (r = 2, Q = 3) | OK |
| §7.2 | 36 → 323 on 6×6; 40 → 649 at n = 40 | max_component 323, 649 | OK |
| §9 | Python 3.13.18, NumPy 2.5.3, gmpy2 2.3.2, 2 CPUs | MANIFEST environment | OK |
| Thm 3.1 | "smallest instance n = 3, c = 1" | n = 2, c = 1: s = 4, only index 0 evaluated, x₁∧x₂ answered UNSAT | nit |
| §1 | 719 manuscripts, 42 % Lean | external | unverified |

## Macros recomputed

Own script on the raw JSON; ratios from raw times, 8×8 total from the histogram.

| Macro | Mine | macros.tex |
|---|---|---|
| fibreDictRatioMin/Max | 0.02 / 1.23 | 0.02 / 1.23 |
| fibreBucketRatioMin/Max; fibreSubcritMin/Max | 0.38 / 103.6; 0.4 / 103.6 | same |
| fibreDensityMin/Max; fibreInstances; fibreExactChecks | 4 / 24; 19; 3 | same |
| toyRunsCount / toyPointsCertified | 12 / 116 | 12 / 116 |
| directFourteen Sums / ClauseEvals / Minutes / Count | 580 / 18,104,320 / 9.2 / 3108 | same |
| auxTwelve Hone / Htwo / Count / D | −12561/32 / −9664 / 736 / 16 | same |
| compTwoThirds | 83–109 (100); 196–229 (252); 685–1656 (635) | same |
| gridEightTotal / gridEightSizes | 660,647,962,955 / 33 | same |
| histOverPlainMin/Max | 1.2 / 1.3 | 1.2 / 1.3 |
| encRatioMin/Max; encRatioGridMax / RandomMax | 0.6 / 1,083; 134 / 1,083 | same |
| encFortyTwenty Seconds / Nodes / Comp / Ratio | 356 / 97,698 / 649 / 1,083 | same |
| histFortyMs / Entries; gridSixEncComp / Ratio; gridSixSixMs | 329 / 26; 323 / 134; 37.9 | same |
| multirow PB / Two / PBentries / TwoEntries | 43.0 / 50.8 / 37 / 100 | same |
| indepPassed/indepChecks; datasetFiles | 74/74; 64 | same |

45 of 45 agree; the 13 others (aux example and environment strings, gridEightMs) are covered by the byte-identical regeneration of `macros.tex`.

## Specific claims

- **Grid totals = OEIS A006506.** 1234, 55447, 5598861, 1280128950, 660647962955 equal my own bitmask row-DP (full histograms equal too) and the OEIS terms a(4)–a(8). **Verified.**
- **Offspring mean.** 2·(1/r)·2·(1/r) = 4/r². **Correct** (heuristic, as stated).
- **Appendix B.** R/Q ≈ 512ν/√λ; log₂(2Q/R) = ½log₂λ − log₂ν − 8, leading ν/2; T₁ ≈ 256ν²/√λ; T₃ = 1/r ≈ ν²/√λ; T₂ ≈ (ν²/λ^{3/2})(λ/64ν)(ν/2) = ν²/(128√λ); hence 256 : 1 : 1/128 and 1 − 1/257. **Leading-order algebra verified.** Numerically: at **n = 10^100** (λ = 333, ν = 9) the integer choices give **r = 0, Q = 0**; at **n = 10^1000** (λ = 3322, ν = 12) **r = 0**, Q = 2, R = 232 — T₂ and T₃ are undefined. With floors removed, T₁ = −1763 and −611 (2Q < R), T₃ = 4.4 and 2.5, T₂ ≈ 0.01. r ≥ 1 first at λ = 2¹⁶; T₁ > 0 first near λ = 2²⁶ (n ≈ 2^{2^26}). At λ = 2^k, k = 26, 32, 64, 128, 256, 512: T₁/T₃ = 5.8, 48, 144, 196, 224, 239; T₃/T₂ = 464, 341, 205, 164, 146, 138; T₃'s share of the exponent 14.6 %, 2.0 %, 0.69 %, 0.51 %, 0.44 %, 0.42 % (limit 0.39 %). The appendix's caveat that the dropped terms are not small at any writable n is accurate; "about one part in 256" (§5.2, abstract) holds only in the limit.
- **§4.3 runs/points.** 12 runs, 116 points, every point `ok`, 7 runs with H(0) = #F. **Correct for the dataset**; Table 1 disagrees (below).
- **E4 flags.** `brute_force_ok` True for 4×4 (n = 16) and random n = 20, the only n ≤ 20 instances (also reproduced by my own NumPy brute force); `equals_histogram` True for all 15 encodings; `plain_equals_total` and `total_equals_oeis` True throughout. **Verified.**

## Inconsistencies found

1. **Table 1/2 vs text.** §4.3 says Table 1 "lists every run of the toy in the dataset" and totals 12 runs/116 points; `make_tables.py` hard-codes an `order` list that omits `toy_direct_n8_q1` (n = 8, Q = 11, 187 sums, H(0) = 85 = #F), so both tables have 11 rows (105 points).
2. **"packed over dict is 0.02–1.23×" (§6.2).** The macro and the table column are dict/packed; packed/dict would be 0.81–57×. And "does not beat a sparse convolution on any instance" (abstract, §1): packed *is* faster on 5 of 7 real-design instances (dict/packed 1.01–1.23). Honest form: "never by more than 1.23×".
3. **"win by 0.4–103.6×" (§6.2).** bucket/packed; 0.4 (r = 4, 10 rows, k = 20) is a loss for packed, while dict wins 21.6× there. "At most a factor of a few at r = 2": actual max 7.7× (k = 18).
4. **"costs 0.6–1,083× more" (abstract, §1).** 0.6× and 0.8× (b = 1 on 6×6, 5×5) are cheaper than the whole histogram.
5. **§4.3 "q = M ≥ n"** at n = 6: q = 7, M = 15.
6. No untraceable number found apart from the external facts above.

## Reproducibility notes

Audit machine = generation machine (Linux-6.18.44-fc-v180, x86_64, 2 vCPU, 7 GB, Python 3.13.18, NumPy 2.5.3, gmpy2 2.3.2): a clean-directory re-execution, not a cross-platform one. `verify_all.py` 6 m 14 s; `make_tables.py` on a scratch copy: 12 tables + 3 figures identical. `paper/main.pdf` carries the same numbers as `main.tex` and is md5-identical to `Measurement-Indexed-Counting-draft1.pdf`. Wall-clock macros are single runs, not re-timed; their deterministic companions (nodes, sums, lanes, entries) all reproduced.

---

## Author's response to the audit (28 October 2026, after the first draft)

All six prose defects were corrected in the released version of the paper; the data, dataset and code results were not changed.

1. **Table 1/2 omitted `toy_direct_n8_q1`.** `make_tables.py` now lists all 12 runs (116 design points); the text and the tables agree.
2. **Inverted label "packed over dict".** The text now states the ratio as dictionary time over packed time (0.02–1.23) and reads it correctly: the packed form is at best 1.23× faster (on 5 small real-design groups) and at worst 57× slower. The abstract and the contribution list say "at best 1.23× faster and at worst 57× slower" instead of "does not beat on any instance".
3. **"Component methods win by 0.4–103.6×".** Replaced by the better of the two component methods against enumeration (macros `fibreBestSubcritMin/Max`, 14.6–299.9× for k ≥ 20 at r ≥ 3) and "at most 7.7× at r = 2" (macro `fibreBestRtwoMax`).
4. **"costs 0.6–1,083× more".** Now "costs 0.6–1,083× the time of the whole histogram, more than the histogram everywhere except at b = 1".
5. **"q = M ≥ n" at n = 6.** Now "q = ⌊M/Q⌋ = 7 ≥ n with M = 15" (macros from the dataset).
6. **Theorem 3.1's "smallest instance".** Now states the n = 2 case (step 4, only index 0 evaluated, x₁ ∧ x₂ answered UNSAT) alongside the n = 3 running example; the companion note was corrected in the same way.

**Appendix B.** The leading-order ratio is now explicitly "the limit"; the appendix reports the exact terms at λ = 2^k (T₁/T₃ = 6, 48, 144, 196, 239 and T₃'s share 14.6 %, 2.0 %, 0.69 %, 0.51 %, 0.42 % for k = 26, 32, 64, 128, 512; limit 0.39 %), that r = 0 at n = 10^100 and 10^1000 under the integer choices, and that T₁ > 0 first at λ = 2^26. These numbers are generated by `make_tables.py` from the parameter formulas and were recomputed independently by the author (same values as the auditor's).

**Process points.** `independent_check.py` no longer writes into the released dataset unless `--write` is passed (it writes `independent_check_rerun.json` in the working directory), and `verify_all.py` no longer exempts any file from the manifest check. Section 9 of the paper now states which runs quick mode skips and what each rests on. The README's runtime estimate was corrected to the measured ~7 minutes.

A second pass by the same auditor on the corrected text follows.

---

## Second pass (auditor)

Scope: the corrected text only against the unchanged dataset (all 64 manifest hashes still match; `independent_check.json` sha `9aaa2bba…`, untouched).

1. **Tables 1–2.** Now 12 rows; the "(2.6) certified" column sums to 116 = `\toyPointsCertified`; the new row (n = 8, Q = 11, 187 sums, 85 = #F, 3.6 s) matches `toy_direct_n8_q1`. §4.3 "lists every run" is now true.
2. **Corrected sentences vs data** (recomputed from raw JSON):
   - dict/packed 0.02–1.23 read as "at best 1.23× faster (5 real-design groups), at worst 57× slower": max packed/dict = 1/0.01743 = 57.4 → 57; exactly 5 instances have dict slower than packed (real design, k = 12–20). Correct.
   - "better component method vs enumeration 14.6–299.9× for k ≥ 20 at r ≥ 3": bucket/min(dict, packed) on the 8 synthetic instances = 102.2, 95.9, 103.6, 58.2, 299.9, 21.6, 14.6, 213.2; "at most 7.7× at r = 2" (k = 18). Correct.
   - "0.6–1,083× the time of the whole histogram, more than the histogram everywhere except at b = 1": every b ≥ 2 encoding is ≥ 4.4×; b = 1 gives 1.1, 0.8, 0.6. Consistent.
   - §4.3 q = 7, M = 15, 1,458 sums; Theorem 3.1 n = 2: step 4, index {0}; n = 3: step 5, indices {0, 5}; §1/§5.2 shares 14.6 %, 2.0 %, 0.69 %. Correct.
3. **Appendix B macros.** All 20 new macros (`expFirstPositive` 26; `expRatioA–E` 6, 48, 144, 196, 239; `expShareA–E` 14.6, 2.0, 0.69, 0.51, 0.42; limit 0.39; `auxSix*`, `fibre*`) equal my own evaluation with exact integer parameters. One wording imprecision remains: "T₁ first becomes positive at λ = 2²⁶" (Appendix B, author's response) — the exact threshold is λ = 44,318,976 ≈ 2^25.4, where T₁ ≈ 0 and T₃'s share is ≈ 100 %; §5.2's "smallest power of two" is the right phrasing, and §1 bullet 2's "smallest n for which its parameters are all defined and positive" is a different condition (r, Q ≥ 1 already at λ = 2¹⁶, where T₁ = −128). The 14.6 % figure is correctly tied to λ = 2²⁶ wherever λ is stated.
4. **Re-execution.** `python3 code/verify_all.py`: `manifest: 64 files, 0 problems`; `regeneration: 0 differences in deterministic quantities`; `74 of 74 independent checks passed`; `VERIFY_ALL: PASS`; 6 m 56 s. The independent check now writes `independent_check_rerun.json` to a temporary cwd and leaves the dataset alone. `make_tables.py` regenerates all 12 table files byte-identically (78 macros); `main.pdf` carries the corrected numbers.

Verdict: **PASS**.

---

**After the second pass (author).** The one remaining imprecision was corrected: Appendix B now gives the exact threshold at which $T_1$ becomes positive ($\lambda = 44{,}318{,}976 \approx 2^{25.4}$, generated as macros `expThresholdLambda`/`expThresholdLog` by bisection on the integer parameter formulas) and calls $\lambda = 2^{26}$ the smallest power of two with $T_1 > 0$; the contribution list ties the 14.6 % figure to $\lambda = 2^{26}$ with that wording. No data changed; `make_tables.py` regenerated (80 macros) and the paper recompiled.
