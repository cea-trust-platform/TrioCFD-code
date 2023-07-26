#!/bin/bash

DOC="plan-de-dev-TrioCFD_2020_2025"

# Generation of the pdf report
pdflatex ${DOC}.tex
bibtex ${DOC}
pdflatex ${DOC}.tex
pdflatex ${DOC}.tex

# Installation of the final pdf report
cp ${DOC}.pdf ../../doc

# Cleaning
for ext in aux bbl blg out log toc pdf
do
  rm "${DOC}.${ext}"
done

rm annexe.aux

exit 0
