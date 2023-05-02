#!/bin/bash

DOC="TrioCFD_PGC"

# Generation of the pdf report
pdflatex ${DOC}.tex
bibtex ${DOC}
pdflatex ${DOC}.tex
pdflatex ${DOC}.tex

# Installation of the final pdf report
cp ${DOC}.pdf ../../doc

# Cleaning
for ext in aux bbl blg idx lof log lot toc pdf
do
  rm "${DOC}.${ext}"
done

exit 0
