#!/bin/bash

DOC="TrioCFD_PGC"

# Generation of the pdf report
pdflatex ${DOC}.tex
bibtex ${DOC}
pdflatex ${DOC}.tex
pdflatex ${DOC}.tex

# Installation of the final pdf report
mv ${DOC}.pdf ../../doc

exit 0
