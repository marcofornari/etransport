```@meta
CurrentModule = Mstar2t
DocTestSetup = quote
    using Mstar2t
end
```

# Documentation

## Types
```@docs
ParabBand(mstar,ϵ₀,type,deg)
```
```@docs
BandStructure(n,bands,μ)
```

## Functions

```@docs
electrical_conductivity(bandstructure::BandStructure, temperature, scattering_time; exportasdf::Bool, fulltensor::Bool)
```
```@docs
seebeck_coefficient(bandstructure::BandStructure, temperature, scattering_time; exportasdf::Bool, fulltensor::Bool)
```
```@docs
thermal_conductivity(bandstructure::BandStructure, temperature, scattering_time; exportasdf::Bool, fulltensor::Bool)
```
```@docs
carrier_concentration(bandstructure::BandStructure, temperature, fermi_level; exportasdf::Bool)
```
```@docs
lorenz_tensor(bandstructure::BandStructure, temperature, scattering_time; exportasdf::Bool, fulltensor::Bool)
```
```@docs
Ln(ϵ₀, μ, bandtype, T, β, n, idx, τ)
```
```@docs
constant(A=1.0)
```
```@docs
T_fun(τ::Function)
```
```@docs
impurity(ϵ_im::Real, A_im::Real=1;γ::Real=1.)
```
```@docs
acoustic(bands_min, A_sm, τm_max; T₀, μ_min, μ_max)
```
```@docs
matthiessen(τ_models::Array{ScModel};γ::Float64=-1.)
```
```@docs
plot(num_plots, x, y, z=Float64[]; titles, xlabels, ylabels, zlabel, color, colorscheme, vlines, annotations)

```@docs
plot_bandstructure(bs,xaxis=range(-1,1,length=100); colors=nothing)
```
```@docs
plot_τ(scm, t, type="e-T"; μ, ϵ₀, bandtype) 
```
```@docs
savedata(path, data)
```
```@docs
savefig(fullpath, fig)
```

## Index

```@index
```