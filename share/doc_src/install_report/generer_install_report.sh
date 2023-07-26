#!/bin/bash

# Generation install_report
DOC="TrioCFD_install_report"

# Generation of the pdf report
pdflatex ${DOC}.tex
pdflatex ${DOC}.tex
pdflatex ${DOC}.tex

# Installation of the final pdf report
#cp ${DOC}.pdf ../../doc

# Cleaning
for ext in aux log out
do
  rm "${DOC}.${ext}"
done

exit 0
