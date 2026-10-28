#!/bin/sh
# Build the release archive and the arXiv source bundle.
set -e
cd "$(dirname "$0")"
rm -rf release && mkdir -p release/arxiv
cp paper/main.tex paper/refs.bib paper/main.bbl release/arxiv/ 2>/dev/null || cp paper/main.tex paper/refs.bib release/arxiv/
cp -r paper/tables paper/figures release/arxiv/
(cd release/arxiv && zip -qr ../measurement-indexed-counting-arxiv.zip .)
zip -qr release/measurement-indexed-counting-repo.zip README.md LICENSE VERIFICATION.md make_release.sh code dataset dataset_meta community paper/main.tex paper/refs.bib paper/main.pdf paper/tables paper/figures -x "paper/tables/.*"
ls -la release
