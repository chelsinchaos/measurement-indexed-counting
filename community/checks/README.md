# checks/ — executable model of the counting pipeline (family 378)

`measurement_toy.py` (Python 3.10+, standard library only, exact arithmetic) implements the clause polynomials,
replacement linear forms, measurement designs (Definition 2.1, condition (4) verified at every point), identity (2.6),
the restricted sums and identity (2.9), Hermite interpolation, and three design families including the paper's own
Section 2.7 construction at small parameters. Every stage is checked against brute-force enumeration.

```
python3 measurement_toy.py --n 12 --m 12 --q 1              # (2.6) at 17 points, (2.9) with q = 1, H(0) = #F = 736
python3 measurement_toy.py --n 12 --m 12 --design aux       # Section 2.7 design at r = 2, Q = 2, R = 2, printed and certified
python3 measurement_toy.py --n 6 --m 6 --seed 2 --design aux   # full pipeline through the Section 2.7 design, H(0) = #F = 34
python3 measurement_toy.py --n 8 --m 8 --design aux --Q 3 --R 4   # the degenerate R = Q + 1 case (no measurements)
```

The full report, dataset (every instance, design and result, with SHA-256 manifest), the independent re-implementation
check and the audit report are in the repository linked from NOTES.md.
