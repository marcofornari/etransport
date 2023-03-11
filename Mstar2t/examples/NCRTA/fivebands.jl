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
import CairoMakie
using LaTeXStrings
using Mstar2t: Scattering


savepath = joinpath("examples", "NCRTA")
num_plots = 3

# BAND STRUCTURE DEFINITION
m_1 = [.8, .8, .8, 0.0, 0.0, 0.0];
ϵ₀_1 = .5;
type_1 = 1;
band_1 = ParabBand(m_1,ϵ₀_1,type_1,1);   # create the band

m_2 = [.5, .5, .5, 0.0, 0.0, 0.0];
ϵ₀_2 = .35;
type_2 = 1;
band_2 = ParabBand(m_2,ϵ₀_2,type_2,1);   # create the band

m_3 = [1., 1., 1., 0.0, 0.0, 0.0];
ϵ₀_3 = .32;
type_3 = -1;
band_3 = ParabBand(m_3,ϵ₀_3,type_3,1);   # create the band

m_4 = [.5, .5, .5, 0.0, 0.0, 0.0];
ϵ₀_4 = .25;
type_4 = -1;
band_4 = ParabBand(m_4,ϵ₀_4,type_4,1);   # create the band

m_5 = [.7, .7, .7, 0.0, 0.0, 0.0];
ϵ₀_5 = .22;
type_5 = -1;
band_5 = ParabBand(m_5,ϵ₀_5,type_5,1);   # create the band

μ = collect(.2:0.01:.6);
model = BandStructure(5,[band_1,band_2,band_3,band_4,band_5],μ);   # build the band structure

T = collect(150:50:650);

T_labels = [L"$T$ $[K]$" for i in 1:num_plots];
μ_labels = [L"$μ$ $[eV]$" for i in 1:num_plots];

ylabels = [L"$\sigma$ $[(\Omega m)^{-1}]$",L"$S$ $[\mu VK^{-1}]$",L"n"];

## 1. constant
τ_form = Scattering.constant()

σ = electrical_conductivity(model,T,τ_form);
S = seebeck_coefficient(model,T,τ_form);
n = carrier_concentration(model,T,τ_form);

T_titles = [L"$σ$ vs $T$, $τ = const$",L"$S$ vs $T$, $τ = const$",L"$n$ vs $T$, $τ = const$"];
fig_T = plot(3, T, [σ,S*1e6,n], μ, titles=T_titles, xlabels=T_labels, ylabels=ylabels, zlabel=μ_labels[1]);
savefig(joinpath(savepath, "ex_fivebands_const_T.png"), fig_T);
μ_titles = [L"$σ$ vs $T$, $τ = const$",L"$S$ vs $T$, $τ = const$",L"$n$ vs $T$, $τ = const$"];
fig_mu = Mstar2t.plot(num_plots, μ, [σ,S*1e6,n], T, titles=μ_titles, xlabels=μ_labels, ylabels=ylabels, zlabel=T_labels[1]);
savefig(joinpath(savepath, "ex_fivebands_const_mu.png"), fig_mu);

## 2. acoustic
τ_form = Scattering.acoustic(model);
# τ 1c
fig = plot_τ(τ_form,T,"μ-T",μ=μ,ϵ₀=ϵ₀_1,bandtype=type_1);
savefig(joinpath(savepath, "tau_acoustic_1c.png"), fig);
# τ 2c
fig = plot_τ(τ_form,T,"μ-T",μ=μ,ϵ₀=ϵ₀_2,bandtype=type_2);
savefig(joinpath(savepath, "tau_acoustic_2c.png"), fig);
# τ 1v
fig = plot_τ(τ_form,T,"μ-T",μ=μ,ϵ₀=ϵ₀_3,bandtype=type_3);
savefig(joinpath(savepath, "tau_acoustic_1v.png"), fig);
# τ 2v
fig = plot_τ(τ_form,T,"μ-T",μ=μ,ϵ₀=ϵ₀_4,bandtype=type_4);
savefig(joinpath(savepath, "tau_acoustic_2v.png"), fig);
# τ 3v
fig = plot_τ(τ_form,T,"μ-T",μ=μ,ϵ₀=ϵ₀_5,bandtype=type_5);
savefig(joinpath(savepath, "tau_acoustic_3v.png"), fig);

σ = electrical_conductivity(model,T,τ_form);
S = seebeck_coefficient(model,T,τ_form);
n = carrier_concentration(model,T,τ_form);

T_titles = [L"$σ$ vs $T$, $τ = acoustic$",L"$S$ vs $T$, $τ = acoustic$",L"$n$ vs $T$, $τ = acoustic$"];
fig_T = plot(3, T, [σ,S*1e6,n], μ, titles=T_titles, xlabels=T_labels, ylabels=ylabels, zlabel=μ_labels[1]);
savefig(joinpath(savepath, "ex_fivebands_acoustic_T.png"), fig_T);
μ_titles = [L"$σ$ vs $T$, $τ = acoustic$",L"$S$ vs $T$, $τ = acoustic$",L"$n$ vs $T$, $τ = acoustic$"];
fig_mu = Mstar2t.plot(num_plots, μ, [σ,S*1e6,n], T, titles=μ_titles, xlabels=μ_labels, ylabels=ylabels, zlabel=T_labels[1]);
savefig(joinpath(savepath, "ex_fivebands_acoustic_mu.png"), fig_mu);
