#!/bin/bash

DOC="plan-de-dev-TrioCFD_2020_2025"

# Generation of the pdf report
pdflatex ${DOC}.tex
bibtex ${DOC}
pdflatex ${DOC}.tex
pdflatex ${DOC}.tex

# Installation of the final pdf report
mv ${DOC}.pdf ../../doc

exit 0
