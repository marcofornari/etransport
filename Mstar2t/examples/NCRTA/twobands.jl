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
num_plots = 3

# BAND STRUCTURE DEFINITION
m_1 = [.5, .5, .5, 0.0, 0.0, 0.0];
ϵ₀_1 = .3;
type_1 = 1;
deg_1 = 1;
band_1 = ParabBand(m_1,ϵ₀_1,type_1,deg_1)   # create the conduction band

m_2 = [1.0, 1.0, 1.0, 0.0, 0.0, 0.0];
ϵ₀_2 = .0;
type_2 = -1;
deg_2 = 1;
band_2 = ParabBand(m_2,ϵ₀_2,type_2,deg_2)   # create the valence band

μ = collect(.2:0.01:.4);
model = BandStructure(2,[band_1,band_2],μ)   # build the two-band structure

T = collect(300:10:650);

T_labels = [L"$T$ $[K]$" for i in 1:num_plots];
μ_labels = [L"$μ$ $[eV]$" for i in 1:num_plots];

ylabels = [L"$\sigma$ $[(\Omega m)^{-1}]$",
L"$S$ $[\mu VK^{-1}]$",
L"n"];

## 1. constant
τ_form = Scattering.constant();

σ = electrical_conductivity(model,T,τ_form);
S = seebeck_coefficient(model,T,τ_form);
n = carrier_concentration(model,T,τ_form);

titles = [L"$σ$ vs $T$, $τ = const$",L"$S$ vs $T$, $τ = const$",L"$n$ vs $T$, $τ = const$"];
fig = plot(3, T, [σ,S*1e6,n], μ, titles=titles, xlabels=T_labels, ylabels=ylabels, zlabel=μ_labels[1]);
savefig(joinpath(savepath, "ex_two_isobands_constant.png"), fig);

## 2. acoustic
τ_form = Scattering.acoustic(model,T₀=180,μ_min=5,μ_max=5);
# conduction band
fig = plot_τ(τ_form,T,"μ-T", μ=μ, ϵ₀=ϵ₀_1, bandtype=type_1);
savefig(joinpath(savepath, "acoustic_tau_cond.png"), fig);
# valence band
fig = plot_τ(τ_form,T,"μ-T", μ=μ, ϵ₀=ϵ₀_2, bandtype=type_2);
savefig(joinpath(savepath, "acoustic_tau_val.png"), fig);

σ_nc = electrical_conductivity(model,T,τ_form);
S_nc = seebeck_coefficient(model,T,τ_form);
n_nc = carrier_concentration(model,T,τ_form);

T_titles = [L"$σ$ vs $T$, $τ = acoustic$",L"$S$ vs $T$, $τ = acoustic$",L"$n$ vs $T$, $τ = acoustic$"];
fig_T = plot(3, T, [σ,S*1e6,n], μ, titles=T_titles, xlabels=T_labels, ylabels=ylabels, zlabel=μ_labels[1]);
savefig(joinpath(savepath, "ex_twobands_acoustic_T.png"), fig_T);
μ_titles = [L"$σ$ vs $T$, $τ = acoustic$",L"$S$ vs $T$, $τ = acoustic$",L"$n$ vs $T$, $τ = acoustic$"];
fig_mu = Mstar2t.plot(num_plots, μ, [σ,S*1e6,n], T, titles=μ_titles, xlabels=μ_labels, ylabels=ylabels, zlabel=T_labels[1]);
savefig(joinpath(savepath, "ex_twobands_acoustic_mu.png"), fig_mu);