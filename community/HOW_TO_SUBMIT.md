# How to publish the three pieces

## 1. arXiv (the report)

- Category: cs.CC (Computational Complexity); cross-list cs.DS and cs.LO.
- Upload `release/measurement-indexed-counting-arxiv.zip` (main.tex, refs.bib, main.bbl, tables/, figures/). arXiv compiles
  with pdflatex; the bundle includes the .bbl so no BibTeX run is needed.
- Title: Measurement-Indexed Counting: An Executable Study of OpenAI's #3-SAT Algorithm, Semiring Component Caches,
  and the Limits of Oblivious Sampling
- Comments field: "17 pages, 3 figures. Code, dataset and verification logs: <repository URL>."
- Abstract: paste from paper/main.tex (the \begin{abstract} block).
- A first submission to cs.CC may require an endorsement; arXiv tells you at upload time.
- Before uploading, replace "<repository URL>" everywhere (README.md, NOTES.md, forum_post.md, the paper's Section 9)
  with the public URL of the repository.

## 2. The community fork (mathcommunity/oai-math-revised)

1. Fork https://github.com/mathcommunity/oai-math-revised to your account (a pseudonymous account is allowed by its README).
2. In your fork, under `preprints/The-Exponential-Time-Hypothesis-is-False-October-25-2026/`, add
   `NOTES.md` and the folder `checks/` (measurement_toy.py, README.md) from `community/`.
3. Open a pull request with the text of `community/PR_DESCRIPTION.md`. The README says all submissions are reviewed.

## 3. The OpenAI Developer Community thread

Thread: "First look at mathematics manuscripts from an internal frontier model at OpenAI"
(https://community.openai.com/t/first-look-at-mathematics-manuscripts-from-an-internal-frontier-model-at-openai/1403886).
Post `community/forum_post.md` as a reply, with the repository URL filled in. Other replies there are verifications of
the same kind (Lean re-runs, numerical checks), so the register fits.

## What not to claim, anywhere

No verification of the theorem beyond toy scale; no speed-up over production counters; no FHC result; no P vs NP
statement. The paper's "What this paper does not claim" paragraph is the reference.
