# Mstar2t Julia package

*Julia package version of the m\*2T software*

## Requirements 

- Julia>=1.6

## Dependencies

- ArgParse == 1.1
- CSV == 0.10
- CairoMakie == 0.10
- Colors == 0.4
- DataFrames == 1.4
- Distributions == 0.25
- Einsum == 0.4
- FastGaussQuadrature == 0.5
- GLMakie == 0.8
- HTTP == 1.7
- HypergeometricFunctions == 0.3
- JSON == 0.21
- LaTeXStrings == 1.3
- LinearAlgebra == 3.4
- Parameters == 0.12
- PlotUtils == 1.3
- PolyLog == 2.3
- Polynomials == 3.2
- QuadGK == 2.6
- SpecialFunctions == 2.1

## Installation

Open a Julia REPL, enter the Julia package manager (typing `]`) and then add the package with the ssh url of the repository:

```bash
(v1.6) pkg> add https://github.com/marcofornari/etrasport.git
```

Then you should be able to import Mstar2t as a normal Julia package and use its functions.  

**Parallel version**: to run the code in parallel with Julia's multi-threading, launch a REPL specifying the number `N` of execution threads you need:
```bash
$ julia -t N

julia> Threads.nthreads() # to check threads available
N
```
Note: By running `julia -t auto`, Julia will use the number of local CPU threads.

### Examples
Check [examples](https://github.com/marcofornari/etrasport/tree/main/julia/examples) folder for examples. 
