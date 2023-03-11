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
using LaTeXStrings
using Mstar2t: Scattering

# BAND STRUCTURE DEFINITION
m_1 = [0.1, 0.1, 0.1, 0.0, 0.0, 0.0];
ϵ₀_1 = 1.0;
type_1 = 1;  # conduction band
deg_1 = 1;
band_1 = ParabBand(m_1,ϵ₀_1,type_1,deg_1);  # create the conduction band

m_2 = [1.0, 0.5, 0.5, 0.0, 0.0, 0.0];
ϵ₀_2 = collect(0.0:0.005:1.05);
type_2 = -1;    # valence band
deg_2 = 1;
band_2 = ParabBand(m_2,ϵ₀_2,type_2,deg_2);  # create the valence band

μ = 0.8;
model = BandStructure(2,[band_1,band_2],μ); # build the two-band structure

T = collect(50:50:700);

τ_form = Scattering.constant();             # relaxation time

# TENSORS COMPUTATION
σ = electrical_conductivity(model,T,τ_form);
S = seebeck_coefficient(model,T,τ_form);
n = carrier_concentration(model,T,τ_form);

# PLOT
gap = ϵ₀_1 .- ϵ₀_2;

titles = [L"$σ$ vs $T$, $μ$ = %$μ, $τ = const$",
L"$S$ vs $T$, $μ$ = %$μ, $τ = const$",
L"$n$ vs $T$, $μ$ = %$μ, $τ = const$"];

xlabels = [L"$T$ $[K]$", 
L"$T$ $[K]$",
L"$T$ $[K]$"];

ylabels = [L"$\sigma$ $[(\Omega m)^{-1}]$",
L"$S$ $[\mu VK^{-1}]$",
L"n"];

zlabel = L"$Gap$ $[eV]$";

savepath = joinpath("examples", "CRTA", "BipolarTransport");

fig = plot(3, T, [σ,S*1e6,n], gap, titles=titles, xlabels=xlabels, ylabels=ylabels, zlabel=zlabel)

savefig(joinpath(savepath, "ex_two_isobands.png"), fig);


