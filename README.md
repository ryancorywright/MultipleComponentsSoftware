# MultipleComponentsSoftware

Software supplement for the paper

"Sparse PCA with Multiple Principal Components"

by Ryan Cory-Wright and Jean Pauphilet for which a preprint is available at https://optimization-online.org/2022/09/sparse-pca-with-multiple-components/

## Introduction

The software in this package is designed to provide certifiably near-optimal solutions to the problem

`max tr(U'QU`
`s.t. ||U||_0 <=k, U'U=I`

using a relax-and-round, exact, or Lagrangean relaxation approach.


## Installation and set up

In order to run this software, you must install a recent version of Julia from http://julialang.org/downloads/, and a recent version of the Mosek solver (academic licenses are freely available at https://www.mosek.com/products/academic-licenses/). To run parts of the code you will also need to install the Gurobi solver (academic licenses are freely available at https://www.gurobi.com/academia/academic-program-and-licenses/) If you do not have access to Mosek, you could use the freely available SCS package instead, although your results may vary.  The most recent version of Julia at the time this code was last tested was Julia 1.7.2 using Mosek (beta) version 10.

Several packages must be installed in Julia before the code can be run.  These packages can be found in "core_julia1p7.jl"

At this point, the "createTablex.jl" and "createFigurex.jl" files should run successfully.  To run the script, navigate to the folder directory and run:

`include("createTablex.jl")`

 The script will reproduce table x in the paper, where x is an integer between 1 and 4.


 Note that the createFigure2.jl script also runs some code originally written by Berk and Bertsimas for benchmarking purposes; to fully replicate our results however, you will also need to install the elasticnet R package and download Matlab to run the code released by Hein and Buhler.


## Citing MultiplePCs.jl

If you use MultiplePCs.jl, we ask that you please cite the following preprint https://optimization-online.org/2022/09/sparse-pca-with-multiple-components/

## Thank you

Thank you for your interest in MultiplePCs. Please let us know if you encounter any issues using this code, or have comments or questions.  Feel free to email us anytime.


Ryan Cory-Wright
r.cory-wright@imperial.ac.uk

Jean Pauphilet
jpauphilet@london.edu
