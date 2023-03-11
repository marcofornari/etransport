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
import CairoMakie

# BAND STRUCTURE DEFINITION
m_1 = [1.0, 1.0, 1.0, 0.0, 0.0, 0.0];
ϵ₀_1 = 0.5;
type_1 = 1;
deg_1 = 1;
band_1 = ParabBand(m_1,ϵ₀_1,type_1,deg_1);   # create the conduction band

m_2 = [1.0, 1.0, 1.0, 0.0, 0.0, 0.0];
ϵ₀_2 = 0.0;
type_2 = -1;
deg_2 = 1;
band_2 = ParabBand(m_2,ϵ₀_2,type_2,deg_2);   # create the "fixed" valence band

m_3 = [1.0, 1.0, 1.0, 0.0, 0.0, 0.0];
ϵ₀_3 = collect(-0.5:0.001:0.5);
type_3 = -1;
deg_3 = 1;
band_3 = ParabBand(m_3,ϵ₀_3,type_3,deg_3);   # create the "moving" valence band

μ = -0.1;
model = BandStructure(3,[band_1,band_2,band_3],μ);   # build the three-band structure

T = collect(50.:50:700);

τ_form = Scattering.constant()

# TENSORS COMPUTATION
σ = electrical_conductivity(model,T,τ_form);
S = seebeck_coefficient(model,T,τ_form);
n = carrier_concentration(model,T,τ_form);

# PLOT
gap = ϵ₀_2 .- ϵ₀_3;

titles = [L"$σ$ vs gap, $μ$ = %$μ, $τ = const$",
          L"$S$ vs gap, $μ$ = %$μ, $τ = const$",
          L"$n$ vs gap, $μ$ = %$μ, $τ = const$"];

xlabels = [L"$Gap$ $[eV]$" for i in 1:3];

ylabels = [L"$\sigma$ $[(\Omega m)^{-1}]$",
           L"$S$ $[\mu VK^{-1}]$",
           L"n"];

zlabel = L"$T$ $[K]$";

savepath = joinpath("examples", "CRTA", "BandConvergence");

# colormap
mypalette = CairoMakie.cgrad(:viridis, 6, categorical=true);
# annotations
annotations = Dict(L"μ" => [(0.11,6e5),(0.11,0),(0.11,5e6)], 
                   L"E_V" => [(0.01,6e5-5e3),(0.01,0),(0.01,5e6-5e4)]);

# transpose the tensors to have the temperature as the z axis and the band gap as the x axis
σ = convert(Matrix{Float64}, transpose(σ));
S = convert(Matrix{Float64}, transpose(S));
n = convert(Matrix{Float64}, transpose(n));
                   
fig = plot(3, gap, [σ,S*1e6,n], T, titles=titles, xlabels=xlabels, ylabels=ylabels, zlabel=zlabel, vlines=[0.0,0.1], color=mypalette[3], annotations=annotations);

savefig(joinpath(savepath, "ex_three_isobands.png"), fig);

