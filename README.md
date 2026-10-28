# Measurement-Indexed Counting

Repository: https://github.com/chelsinchaos/measurement-indexed-counting

Code, dataset and paper for *Measurement-Indexed Counting: An Executable Study of OpenAI's #3-SAT Algorithm, Semiring Component Caches, and the Limits of Oblivious Sampling* (C. A. McElveen, October 2026).

```
code/      measurement_toy.py    exact executable model of the OpenAI #3-SAT pipeline (stdlib only)
           fibre_loop.py         the fibre loop of identity (2.6) three ways (numpy, gmpy2)
           constrained_count.py  histogram-valued component cache vs cardinality encodings (stdlib only)
           run_all.py            regenerates the dataset (verification run 1)
           independent_check.py  recomputes every result by different algorithms (verification run 2)
           verify_all.py         one command: regenerate (--quick), compare, independent check, manifest
           make_tables.py        produces every table, figure and number of the paper from the dataset
dataset/   instances/ (DIMACS), designs/ (JSON), results/ (JSON), logs/, MANIFEST.json (SHA-256)
paper/     main.tex, refs.bib, tables/ and figures/ (generated), main.pdf
community/ NOTES.md, forum_post.md, PR_DESCRIPTION.md for the community revision of family 378
VERIFICATION.md   report of the independent audit (verification run 3)
```

## Reproduce

```
cd code
python3 verify_all.py                 # ~7 min: regenerate (quick), compare, independent check, manifest
python3 run_all.py                    # ~30 min: full regeneration including the two long runs
python3 make_tables.py && cd ../paper && latexmk -pdf main.tex
```

Requirements: Python 3.10+; `numpy` and `gmpy2` for `fibre_loop.py` only; LaTeX for the paper.

## Quick looks

```
python3 code/measurement_toy.py --n 12 --m 12 --q 1              # (2.6), (2.9), H(0) = #F at n = 12
python3 code/measurement_toy.py --n 12 --m 12 --design aux       # the paper's Section 2.7 design, printed
python3 code/constrained_count.py --grids 4x4,5x5,6x6            # grid independent sets by size vs encodings
python3 code/fibre_loop.py --mode synthetic --r 4 --rows 8 --sizes 16,18,20
```

## Licences

Code: MIT (`LICENSE`). Dataset and paper text: CC BY 4.0 (`dataset/LICENSE`). The OpenAI manuscript studied here is distributed by OpenAI under Apache-2.0 at https://github.com/openai/math (family 378).
