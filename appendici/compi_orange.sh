#!/bin/bash

pdflatex $1
makeindex $1.idx -s StyleInd.ist
biber $1
pdflatex $1
pdflatex $1
