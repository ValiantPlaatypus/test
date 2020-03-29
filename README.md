# FluidMechanics-ITA

This repo contains some personal material (notes and some programs) prepared during the teaching activity for the course of Fluid Mechanics, a 3rd-year BSc course in Politecnico di Milano.

So far, the notes are available in Italian only.

## Compiling .tex files
Each of the three main folders contains a compi*.sh file. These files are bash scripts (they must be exectuable!) that drive the compilation of the .pdf documents. Each of these .sh files expect the main .tex file (the one that contains \begin{document} ...) as the only argument, for the compilation of the documents on a Linux platform. As an example, running

 ./compi_orange mainEse
 
in template directory should produce the file mainEse.pdf. The main *.tex files in ./qa and ./appendici are qa.tex and appendici.tex, respectively.

For users relying on a Latex editor (Texmaker, ...), importing main *.tex and running the compilation on this file should be enough.

## Contents
### Notes
This repo contains three main documents: 
- "Notes about Fluid Mechanics";
- "Review of Mathematical Analysis and Introduction to Tensor Algebra and Calculus";
- "Questions and Answers".

### Programs
This repo contains some programs and scripts about:
- Archimedes' principle and buoyant force: equilibrium and stability
- Flow around a 2D airfoil, with panel method: panel method; integration of integral boundary layer equations; force evaluation with control volume methods; uncertainty estimation of the computed loads.
- Integration of the Lorenz's dynamical system, as a reduced model of the Boussinesq equations for natural convection.

## Contributions
Contributions of different nature would be really appreciated, as
- Translation in English of the notes (or part of them)
- Correction of the existing material
- Additional material