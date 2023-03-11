# Usage

There are two options to run Mstar2t: GUI or Command Line Interface (CLI). Both ways are designed to automatically set up the correct environment for the computing server. 

## Usage (with GUI)

```bash
(Interface) $ cd Interface/GUI
(Interface) $ python run_gui.py
```

## Usage (with CLI)

First run the server (computing unit)

```bash
(Interface) $ cd Interface/CLI
(Interface) $ python run_cli.py
```

On a new shell, send simulations requests to the server by running:

```bash
(Interface) $ cd Interface/CLI
(Interface) $ python compute.py -i <input_file> --<tensor_name> --<plot>
```

**Note:** before running a calculation, edit the `results fullpath` argument in the input_file. This path identifies the location where the results are exported and must be in the **same machine** in which the server is running.

## Help (CLI Python interface)

```bash
(Interface) $ python compute.py --help

usage: compute.py [-h] -i INPUTFILE [--conductivity] [--seebeck] [--thermal] [--concentration] [--tplot] [--muplot]

optional arguments:
  -h, --help            show this help message and exit
  -i INPUTFILE, --inputfile INPUTFILE path to input file
  --conductivity, -e    compute electrical conductivity
  --seebeck, -s         compute Seebeck coefficient
  --thermal, -k         compute thermal conductivity
  --concentration, -n   compute carrier concentration
  --tplot               Temperature plot of the tensors
  --muplot              Fermi level plot of the results
```

## Troubleshooting

If the environment gets corrupted and when launching the server you get a `LoadError: ArgumentError` exception similar to this:
```bash
ERROR: LoadError: ArgumentError: <PackageName> is required but does not seem to be installed:
 - Run `Pkg.instantiate()` to install all recorded dependencies.
 ```
 you may try the following steps (see below):
 1. Remove the `Manifest.toml` file inside `Mstar2t`;
 2. Open the julia REPL;
 3. Enter the julia package manager (typing `]`); 
 4. Activate the `Mstar2t` env;
 5. Run `resolve`;
 6. Run `instantiate`.

```bash
$ cd Mstar2t
$ rm Manifest.toml 
$ julia               _

_       _ _(_)_     |  Documentation: https://docs.julialang.org
  (_)     | (_) (_)    |
   _ _   _| |_  __ _   |  Type "?" for help, "]?" for Pkg help.
  | | | | | | |/ _` |  |
  | | |_| | | | (_| |  |  Version 1.6.0-rc2 (2021-03-11)
 _/ |\__'_|_|_|\__'_|  |  Official https://julialang.org/ release
|__/                   |

(@v1.6) pkg> activate .
(Mstar2t) pkg> resolve
(Mstar2t) pkg> instantiate
 ```
This should clean and rebuild the environment with the correct dependencies. 
