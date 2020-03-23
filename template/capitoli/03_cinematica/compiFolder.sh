#!/bin/bash

shopt -s nullglob
for fname in *.tex; do 
   if [[ $fname != *in.tex && $fname != "logicNumb.tex" 
           && $fname != "aerodinamica.tex" ]]; then
      echo
      echo " -------------- My name is $fname -----------------------"
      pdflatex -output-directory $(pwd) $(pwd)/$fname
      echo
      echo
   fi
done
rm *~
mv *.aux ./aux
mv *.log ./log
mv *.pdf ./pdf
