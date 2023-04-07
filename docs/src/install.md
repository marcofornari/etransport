# Installation

## Requirements 

- Julia >= 1.6
- python 3

## Dependencies

### Computing unit (Julia)

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

### Interface

- requests
- numpy
- pandas
- matplotlib
- colorama
- argparse
- PySide2==5.15

## Installation

HTTPS: 
```bash	
$ git clone https://github.com/marcofornari/etrasport.git
```

SSH:
```bash	
$ git clone git@github.com:marcofornari/etrasport.git
```

### Computing unit

```bash
$ cd Mstar2t
$ julia
```

Enter the Julia package manager (typing `]`) and run the following commands:

```bash
(v1.6) pkg> activate .
(Mstar2t) pkg> instantiate
```

### Installation (Interface)

Unix/Mac:

```bash
$ cd Interface
$ python3.7 -m venv Interface
$ source Interface/bin/activate
(Interface) $ pip install --upgrade pip
(Interface) $ python -m pip install -r requirements.txt
```

Windows:

```Powershell
> cd Interface
> python3.7 -m venv Interface
> .\Interface\Scripts\activate
(Interface) > pip install --upgrade pip
(Interface) > python -m pip install -r requirements.txt
```
