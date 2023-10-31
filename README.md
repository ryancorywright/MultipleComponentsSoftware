# MultipleComponentsSoftware

Software supplement for the paper

"Sparse PCA with Multiple Principal Components"

by Ryan Cory-Wright and Jean Pauphilet for which a preprint is available [here](https://optimization-online.org/2022/09/sparse-pca-with-multiple-components/)

## Introduction

The software in this package is designed to provide certifiably near-optimal solutions to the problem

`max tr(U'QU)`
`s.t. ||U||_0 <=k, U'U=I`

using a relax-and-round, exact, or Lagrangian relaxation approach.


## Installation and set up

To run this software, you must install a recent version of Julia from [here](http://julialang.org/downloads/), and a recent version of the Mosek solver (academic licenses are freely available [here](https://www.mosek.com/products/academic-licenses/)). To run parts of the code you will also need to install the Gurobi solver (academic licenses are freely available [here](https://www.gurobi.com/academia/academic-program-and-licenses/)). If you do not have access to Mosek, you could use another semidefinite solver, although your results may vary.  The most recent version of Julia at the time this code was last tested was Julia 1.9.1 using Mosek version 10.1.11 and Gurobi version 10.0.0, the code has previously successfully run using Julia version 1.7. 

Several packages must be installed in Julia before the code can be run.  These packages can be found in "core_julia1p7.jl"

At this point, the "createTablex.jl" and "createFigurex.jl" files should run successfully.  To run the script, navigate to the folder directory and run:

`include("createTablex.jl")`

 The script will reproduce table x in the paper, where x is an integer between 1 and 4.


 Note that the createFigure2.jl script also runs some code originally written by Berk and Bertsimas for benchmarking purposes; to fully replicate our results however, you will also need to install the elasticnet R package and download Matlab to run the code released by Hein and Buhler.

## Correspondence between script files and figures/tables in the paper
- The file "createFigure3.jl" produces the output used to generate Figure 3, Table 3, and Figure EC.1 in our manuscript. 

## Related Packages

If you are interested in computing a single sparse PC, you may want to check out our related package https://github.com/ryancorywright/ScalableSPCA.jl which computes a single provably near-optimal sparse PC via branch-and-cut or relax-and-round techniques.


## Citing MultiplePCs.jl

If you use MultiplePCs.jl, we ask that you please consider citing the following [preprint](https://optimization-online.org/2022/09/sparse-pca-with-multiple-components/):
```
@article{cory2022sparse,
  title={Sparse PCA With Multiple Components},
  author={Cory-Wright, Ryan and Pauphilet, Jean},
  journal={arXiv preprint arXiv:2209.14790},
  year={2022}
}
```

## Thank you

Thank you for your interest in MultiplePCs. Please let us know if you encounter any issues using this code, or have comments or questions.  Feel free to email us anytime.


Ryan Cory-Wright
r [dot] cory-wright [at] imperial [dot] ac [dot] uk

Jean Pauphilet
jpauphilet [at] london [dot] edu
