# *************************************************************************** #
# *                                                                         * #
# *         Mstar2t - Central Michigan University University, 2023          * #
# *                                                                         * #
# *************************************************************************** #
#  This file is part of Mstar2t.                                              #                        
#                                                                             #
#  Mstar2t is free software: you can redistribute it and/or modify it under   #
#  the terms of the GNU General Public License as published by the Free       #
#  Software Foundation, either version 3 of the License, or (at your option)  #
#  any later version.                                                         #
#                                                                             #
#  Mstar2t is distributed in the hope that it will be useful, but WITHOUT     #
#  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or      #
#  FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for   #
#  more details.                                                              #
#                                                                             #
#  You should have received a copy of the GNU General Public License along    #
#  with this program. If not, see <http://www.gnu.org/licenses/>.             #
#                                                                             #
# *************************************************************************** #


using Mstar2t
using Mstar2t: Scattering
using LaTeXStrings


savepath = joinpath("examples", "NCRTA");

# BAND STRUCTURE DEFINITION
m = [0.5, 0.5, 0.5, 0.0, 0.0, 0.0];
ϵ₀ = 0.0;
type = -1;
deg = 1;
band_1 = ParabBand(m,ϵ₀,type,deg);   # create the band

μ = collect(-0.1:0.005:0.2);
model = BandStructure(1,band_1,μ);   # build the band structure

T = collect(51:10:650);

xlabels = [L"$T$ $[K]$", 
L"$T$ $[K]$",
L"$T$ $[K]$"];

ylabels = [L"$\sigma$ $[(\Omega m)^{-1}]$",
L"$S$ $[\mu VK^{-1}]$",
L"n"];

zlabel = L"$μ$ $[eV]$";

## 1. constant
τ_form = Scattering.constant()

σ = electrical_conductivity(model,T,τ_form);
S = seebeck_coefficient(model,T,τ_form);
n = carrier_concentration(model,T,τ_form);

# PLOTS
titles = [L"$σ$ vs $T$, $τ = const$",L"$S$ vs $T$, $τ = const$",L"$n$ vs $T$, $τ = const$"];
fig = plot(3, T, [σ,S*1e6,n], μ, titles=titles, xlabels=xlabels, ylabels=ylabels, zlabel=zlabel);
savefig(joinpath(savepath, "ex_isoband_val_constant.png"), fig);

## 2. acoustic
T₀=45;
μ_min=1;
μ_max=1;
τ_form = Scattering.acoustic(model,1e-1,1e-1,T₀=T₀,μ_min=μ_min,μ_max=μ_max);

σ_nc = electrical_conductivity(model,T,τ_form);
S_nc = seebeck_coefficient(model,T,τ_form);
n_nc = carrier_concentration(model,T,τ_form);

titles = [L"$σ$ vs $T$, $τ = acoustic$",L"$S$ vs $T$, $τ = acoustic$",L"$n$ vs $T$, $τ = acoustic$"];
fig = plot(3, T, [σ_nc,S_nc*1e6,n_nc], μ, titles=titles, xlabels=xlabels, ylabels=ylabels, zlabel=zlabel);
savefig(joinpath(savepath, "ex_isoband_val_acoustic.png"), fig);


## 3. impurity
τ_im = 0.1;
A_im = .5;
γ_im = 0.2;
τ_form = Scattering.impurity(τ_im,A_im,γ=γ_im)

σ_nc = electrical_conductivity(model,T,τ_form);
S_nc = seebeck_coefficient(model,T,τ_form);
n_nc = carrier_concentration(model,T,τ_form);

titles = [L"$σ$ vs $T$, $τ = impurity$",L"$S$ vs $T$, $τ = impurity$",L"$n$ vs $T$, $τ = impurity$"];
fig = plot(3, T, [σ_nc,S_nc*1e6,n_nc], μ, titles=titles, xlabels=xlabels, ylabels=ylabels, zlabel=zlabel);
savefig(joinpath(savepath, "ex_isoband_val_impurity.png"), fig);


## 4. Matthiessen
τ_form1 = Scattering.impurity(τ_im,A_im,γ=γ_im)
τ_form2 = Scattering.acoustic(model,T₀=T₀,μ_min=μ_min,μ_max=μ_max);
τ_form = Scattering.matthiessen([τ_form1,τ_form2],γ=-5.)

σ_nc = electrical_conductivity(model,T,τ_form);
S_nc = seebeck_coefficient(model,T,τ_form);
n_nc = carrier_concentration(model,T,τ_form);

titles = [L"$σ$ vs $T$, $τ = matthiessen$",L"$S$ vs $T$, $τ = matthiessen$",L"$n$ vs $T$, $τ = matthiessen$"];
fig = plot(3, T, [σ_nc,S_nc*1e6,n_nc], μ, titles=titles, xlabels=xlabels, ylabels=ylabels, zlabel=zlabel);
savefig(joinpath(savepath, "ex_isoband_val_matthiessen.png"), fig);
