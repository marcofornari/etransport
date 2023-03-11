var documenterSearchIndex = {"docs":
[{"location":"usage/#Usage","page":"Usage","title":"Usage","text":"","category":"section"},{"location":"usage/","page":"Usage","title":"Usage","text":"There are two options to run Mstar2t: GUI or Command Line Interface (CLI). Both ways are designed to automatically set up the correct environment for the computing server. ","category":"page"},{"location":"usage/#Usage-(with-GUI)","page":"Usage","title":"Usage (with GUI)","text":"","category":"section"},{"location":"usage/","page":"Usage","title":"Usage","text":"(Interface) $ cd Interface/GUI\r\n(Interface) $ python run_gui.py","category":"page"},{"location":"usage/#Usage-(with-CLI)","page":"Usage","title":"Usage (with CLI)","text":"","category":"section"},{"location":"usage/","page":"Usage","title":"Usage","text":"First run the server (computing unit)","category":"page"},{"location":"usage/","page":"Usage","title":"Usage","text":"(Interface) $ cd Interface/CLI\r\n(Interface) $ python run_cli.py","category":"page"},{"location":"usage/","page":"Usage","title":"Usage","text":"On a new shell, send simulations requests to the server by running:","category":"page"},{"location":"usage/","page":"Usage","title":"Usage","text":"(Interface) $ cd Interface/CLI\r\n(Interface) $ python compute.py -i <input_file> --<tensor_name> --<plot>","category":"page"},{"location":"usage/","page":"Usage","title":"Usage","text":"Note: before running a calculation, edit the results fullpath argument in the input_file. This path identifies the location where the results are exported and must be in the same machine in which the server is running.","category":"page"},{"location":"usage/#Help-(CLI-Python-interface)","page":"Usage","title":"Help (CLI Python interface)","text":"","category":"section"},{"location":"usage/","page":"Usage","title":"Usage","text":"(Interface) $ python compute.py --help\r\n\r\nusage: compute.py [-h] -i INPUTFILE [--conductivity] [--seebeck] [--thermal] [--concentration] [--tplot] [--muplot]\r\n\r\noptional arguments:\r\n  -h, --help            show this help message and exit\r\n  -i INPUTFILE, --inputfile INPUTFILE path to input file\r\n  --conductivity, -e    compute electrical conductivity\r\n  --seebeck, -s         compute Seebeck coefficient\r\n  --thermal, -k         compute thermal conductivity\r\n  --concentration, -n   compute carrier concentration\r\n  --tplot               Temperature plot of the tensors\r\n  --muplot              Fermi level plot of the results","category":"page"},{"location":"usage/#Troubleshooting","page":"Usage","title":"Troubleshooting","text":"","category":"section"},{"location":"usage/","page":"Usage","title":"Usage","text":"If the environment gets corrupted and when launching the server you get a LoadError: ArgumentError exception similar to this:","category":"page"},{"location":"usage/","page":"Usage","title":"Usage","text":"ERROR: LoadError: ArgumentError: <PackageName> is required but does not seem to be installed:\r\n - Run `Pkg.instantiate()` to install all recorded dependencies.\r\n ```\r\n you may try the following steps (see below):\r\n 1. Remove the `Manifest.toml` file inside `Mstar2t`;\r\n 2. Open the julia REPL;\r\n 3. Enter the julia package manager (typing `]`); \r\n 4. Activate the `Mstar2t` env;\r\n 5. Run `resolve`;\r\n 6. Run `instantiate`.\r\n","category":"page"},{"location":"usage/","page":"Usage","title":"Usage","text":"bash $ cd Mstar2t $ rm Manifest.toml  $ julia               _","category":"page"},{"location":"usage/","page":"Usage","title":"Usage","text":"_       _ ()_     |  Documentation: https://docs.julialang.org   ()     | () ()    |    _ _   _| |  __ _   |  Type \"?\" for help, \"]?\" for Pkg help.   | | | | | | |/ ` |  |   | | || | | | (| |  |  Version 1.6.0-rc2 (2021-03-11)  _/ |__'|||__'_|  |  Official https://julialang.org/ release |__/                   |","category":"page"},{"location":"usage/","page":"Usage","title":"Usage","text":"(@v1.6) pkg> activate . (Mstar2t) pkg> resolve (Mstar2t) pkg> instantiate  ``` This should clean and rebuild the environment with the correct dependencies. ","category":"page"},{"location":"docfun/","page":"Documentation","title":"Documentation","text":"CurrentModule = Mstar2t\r\nDocTestSetup = quote\r\n    using Mstar2t\r\nend","category":"page"},{"location":"docfun/#Documentation","page":"Documentation","title":"Documentation","text":"","category":"section"},{"location":"docfun/#Types","page":"Documentation","title":"Types","text":"","category":"section"},{"location":"docfun/","page":"Documentation","title":"Documentation","text":"ParabBand(mstar,ϵ₀,type,deg)","category":"page"},{"location":"docfun/","page":"Documentation","title":"Documentation","text":"BandStructure(n,bands,μ)","category":"page"},{"location":"docfun/#Functions","page":"Documentation","title":"Functions","text":"","category":"section"},{"location":"docfun/","page":"Documentation","title":"Documentation","text":"electrical_conductivity(bandstructure::BandStructure, temperature, scattering_time; exportasdf::Bool, fulltensor::Bool)","category":"page"},{"location":"docfun/","page":"Documentation","title":"Documentation","text":"seebeck_coefficient(bandstructure::BandStructure, temperature, scattering_time; exportasdf::Bool, fulltensor::Bool)","category":"page"},{"location":"docfun/","page":"Documentation","title":"Documentation","text":"thermal_conductivity(bandstructure::BandStructure, temperature, scattering_time; exportasdf::Bool, fulltensor::Bool)","category":"page"},{"location":"docfun/","page":"Documentation","title":"Documentation","text":"carrier_concentration(bandstructure::BandStructure, temperature, fermi_level; exportasdf::Bool)","category":"page"},{"location":"docfun/","page":"Documentation","title":"Documentation","text":"lorenz_tensor(bandstructure::BandStructure, temperature, scattering_time; exportasdf::Bool, fulltensor::Bool)","category":"page"},{"location":"docfun/","page":"Documentation","title":"Documentation","text":"Ln(ϵ₀, μ, bandtype, T, β, n, idx, τ)","category":"page"},{"location":"docfun/","page":"Documentation","title":"Documentation","text":"constant(A=1.0)","category":"page"},{"location":"docfun/","page":"Documentation","title":"Documentation","text":"T_fun(τ::Function)","category":"page"},{"location":"docfun/#Mstar2t.Scattering.T_fun-Tuple{Function}","page":"Documentation","title":"Mstar2t.Scattering.T_fun","text":"T_fun(τ::Function)\n\nSet the temperature dependence of the relaxation time from a function given by the user.  \n\nExample: f(T) = sqrt(T)  # τ ∝ √T τmodel = Tfun(f)\n\n\n\n\n\n","category":"method"},{"location":"docfun/","page":"Documentation","title":"Documentation","text":"impurity(ϵ_im::Real, A_im::Real=1;γ::Real=1.)","category":"page"},{"location":"docfun/#Mstar2t.Scattering.impurity","page":"Documentation","title":"Mstar2t.Scattering.impurity","text":"impurity(ϵ_im::Real, A_im::Real=1;γ::Real=1.)\n\nInclude impurity scattering in the simulations.\n\nParameters: ϵim: energy of the impurity in eV A: multiplicative constant in front of the functional expression γ: γ-parameter three-parameter Lorentzian function (Ref: https://en.wikipedia.org/wiki/Cauchydistribution)\n\nExample: ϵim = 0.2  # eV τmodel = impurity(ϵ_im,A=1e-1)\n\n\n\n\n\n","category":"function"},{"location":"docfun/","page":"Documentation","title":"Documentation","text":"acoustic(bands_min, A_sm, τm_max; T₀, μ_min, μ_max)","category":"page"},{"location":"docfun/","page":"Documentation","title":"Documentation","text":"matthiessen(τ_models::Array{ScModel};γ::Float64=-1.)","category":"page"},{"location":"docfun/#Mstar2t.Scattering.matthiessen-Tuple{Array{ScModel, N} where N}","page":"Documentation","title":"Mstar2t.Scattering.matthiessen","text":"matthiessen(τ_models::Array{ScModel};γ::Float64=-1.)\n\nThis function applies Matthiessen's rule to sum up different scattering mechanisms. \n\nParameters: τmodels: vector of relaxation time models γ: exponent in the generalized Matthiessen's rule: τ = (\\sumi τi^γ)^(1/γ)$, where each \\taui$ can be a function of T and/or μ. \n\nExample: ϵim = 0.2  # eV τ1 = constant() τ2 = impurity(ϵim,A=1e-1) τ_model =  matthiessen([τ1,τ2])\n\n\n\n\n\n","category":"method"},{"location":"docfun/","page":"Documentation","title":"Documentation","text":"plot(num_plots, x, y, z=Float64[]; titles, xlabels, ylabels, zlabel, color, colorscheme, vlines, annotations)\r\n","category":"page"},{"location":"docfun/","page":"Documentation","title":"Documentation","text":"@docs plot_bandstructure(bs,xaxis=range(-1,1,length=100); colors=nothing)","category":"page"},{"location":"docfun/","page":"Documentation","title":"Documentation","text":"","category":"page"},{"location":"docfun/","page":"Documentation","title":"Documentation","text":"@docs plot_τ(scm, t, type=\"e-T\"; μ, ϵ₀, bandtype) ","category":"page"},{"location":"docfun/","page":"Documentation","title":"Documentation","text":"","category":"page"},{"location":"docfun/","page":"Documentation","title":"Documentation","text":"@docs savedata(path, data)","category":"page"},{"location":"docfun/","page":"Documentation","title":"Documentation","text":"","category":"page"},{"location":"docfun/","page":"Documentation","title":"Documentation","text":"@docs savefig(fullpath, fig)","category":"page"},{"location":"docfun/","page":"Documentation","title":"Documentation","text":"\r\n## Index\r\n","category":"page"},{"location":"docfun/","page":"Documentation","title":"Documentation","text":"@index ```","category":"page"},{"location":"examples/#Examples","page":"Examples","title":"Examples","text":"","category":"section"},{"location":"examples/","page":"Examples","title":"Examples","text":"Check Usage for examples on how to use the CLI and GUI version of m*2T.","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"The examples folder collects examples of the usage of m*2T as a Julia library. ","category":"page"},{"location":"examples/#Constant-relaxation-time-approximation-CRTA","page":"Examples","title":"Constant relaxation time approximation CRTA","text":"","category":"section"},{"location":"examples/#Single-band-model","page":"Examples","title":"Single-band model","text":"","category":"section"},{"location":"examples/","page":"Examples","title":"Examples","text":"Link","category":"page"},{"location":"examples/#Bipolar-transport","page":"Examples","title":"Bipolar transport","text":"","category":"section"},{"location":"examples/","page":"Examples","title":"Examples","text":"Link","category":"page"},{"location":"examples/#Band-convergence","page":"Examples","title":"Band convergence","text":"","category":"section"},{"location":"examples/","page":"Examples","title":"Examples","text":"Link","category":"page"},{"location":"examples/#Wiedemann-Franz-law","page":"Examples","title":"Wiedemann-Franz law","text":"","category":"section"},{"location":"examples/","page":"Examples","title":"Examples","text":"Link","category":"page"},{"location":"examples/#Non-constant-relaxation-time-approximation-NCRTA","page":"Examples","title":"Non-constant relaxation time approximation NCRTA","text":"","category":"section"},{"location":"examples/#Single-band-model-2","page":"Examples","title":"Single-band model","text":"","category":"section"},{"location":"examples/","page":"Examples","title":"Examples","text":"Link","category":"page"},{"location":"examples/#Double-band-model","page":"Examples","title":"Double-band model","text":"","category":"section"},{"location":"examples/","page":"Examples","title":"Examples","text":"Link","category":"page"},{"location":"examples/#Three-band-model","page":"Examples","title":"Three-band model","text":"","category":"section"},{"location":"examples/","page":"Examples","title":"Examples","text":"Link","category":"page"},{"location":"examples/#Four-band-model","page":"Examples","title":"Four-band model","text":"","category":"section"},{"location":"examples/","page":"Examples","title":"Examples","text":"Link","category":"page"},{"location":"examples/#Five-band-model","page":"Examples","title":"Five-band model","text":"","category":"section"},{"location":"examples/","page":"Examples","title":"Examples","text":"Link","category":"page"},{"location":"install/#Installation","page":"Installation","title":"Installation","text":"","category":"section"},{"location":"install/#Requirements","page":"Installation","title":"Requirements","text":"","category":"section"},{"location":"install/","page":"Installation","title":"Installation","text":"Julia >= 1.6\npython 3","category":"page"},{"location":"install/#Dependencies","page":"Installation","title":"Dependencies","text":"","category":"section"},{"location":"install/#Computing-unit-(Julia)","page":"Installation","title":"Computing unit (Julia)","text":"","category":"section"},{"location":"install/","page":"Installation","title":"Installation","text":"ArgParse == 1.1\nCSV == 0.10\nCairoMakie == 0.10\nColors == 0.4\nDataFrames == 1.4\nDistributions == 0.25\nEinsum == 0.4\nFastGaussQuadrature == 0.5\nGLMakie == 0.8\nHTTP == 1.7\nHypergeometricFunctions == 0.3\nJSON == 0.21\nLaTeXStrings == 1.3\nLinearAlgebra == 3.4\nParameters == 0.12\nPlotUtils == 1.3\nPolyLog == 2.3\nPolynomials == 3.2\nQuadGK == 2.6\nSpecialFunctions == 2.1","category":"page"},{"location":"install/#Interface","page":"Installation","title":"Interface","text":"","category":"section"},{"location":"install/","page":"Installation","title":"Installation","text":"numpy==1.21\npandas==1.3\nrequests==2.24\nPySide2==5.15\nmatplotlib==3.4\ncolorama=0.4\nargparse==1.1","category":"page"},{"location":"install/#Installation-2","page":"Installation","title":"Installation","text":"","category":"section"},{"location":"install/","page":"Installation","title":"Installation","text":"HTTPS: ","category":"page"},{"location":"install/","page":"Installation","title":"Installation","text":"$ git clone https://github.com/marcofornari/etrasport.git","category":"page"},{"location":"install/","page":"Installation","title":"Installation","text":"SSH:","category":"page"},{"location":"install/","page":"Installation","title":"Installation","text":"$ git clone git@github.com:marcofornari/etrasport.git","category":"page"},{"location":"install/#Computing-unit","page":"Installation","title":"Computing unit","text":"","category":"section"},{"location":"install/","page":"Installation","title":"Installation","text":"$ cd Mstar2t\r\n$ julia","category":"page"},{"location":"install/","page":"Installation","title":"Installation","text":"Enter the Julia package manager (typing ]) and run the following commands:","category":"page"},{"location":"install/","page":"Installation","title":"Installation","text":"(v1.6) pkg> activate .\r\n(Mstar2t) pkg> instantiate","category":"page"},{"location":"install/#Installation-(Interface)","page":"Installation","title":"Installation (Interface)","text":"","category":"section"},{"location":"install/","page":"Installation","title":"Installation","text":"Unix/Mac:","category":"page"},{"location":"install/","page":"Installation","title":"Installation","text":"$ cd Interface\r\n$ python3.7 -m venv Interface\r\n$ source Interface/bin/activate\r\n(Interface) $ pip install --upgrade pip\r\n(Interface) $ python -m pip install -r requirements.txt","category":"page"},{"location":"install/","page":"Installation","title":"Installation","text":"Windows:","category":"page"},{"location":"install/","page":"Installation","title":"Installation","text":"> cd Interface\r\n> python3.7 -m venv Interface\r\n> .\\Interface\\Scripts\\activate\r\n(Interface) > pip install --upgrade pip\r\n(Interface) > python -m pip install -r requirements.txt","category":"page"},{"location":"#Mstar2t-Package","page":"Home","title":"Mstar2t Package","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Transport properties from multi-valley density of states in the relaxation time approximation.","category":"page"},{"location":"#Software-Features","page":"Home","title":"Software Features","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Compute electrical conductivity, Seebeck coefficient, carrier concentration, thermal conductivity and Lorentz tensor for multi-valley anisotropic density of states;\nData visualization: transport properties as functions of temperature and chemical potential;\nCommand-line interface and Graphical User Interface;\nMulti-threaded calculations of transport coefficients;\nComparison between experimental and simulated data;\nDatabase creation of transport properties of multi-valley density of states.","category":"page"},{"location":"","page":"Home","title":"Home","text":"Some examples of computations using Mstar2t can be found on the Examples page.","category":"page"},{"location":"","page":"Home","title":"Home","text":"See the Index for the complete list of documented functions and types.","category":"page"},{"location":"#Manual-Outline","page":"Home","title":"Manual Outline","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Pages = [\r\n    \"install.md\",\r\n    \"usage.md\",\r\n    \"docfun.md\",\r\n    \"examples.md\",\r\n    ]\r\nDepth = 1","category":"page"},{"location":"#main-index","page":"Home","title":"Functions","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Pages = [\"docfun.md\"]","category":"page"}]
}
