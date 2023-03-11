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

# BAND STRUCTURE DEFINITION
m_1 = [1.0, 1.0, 1.0, 0.0, 0.0, 0.0];
ϵ₀_1 = collect(0.5:-0.005:-0.25);
type_1 = 1;
deg_1 = 1;

m_2 = [1.0, 1.0, 1.0, 0.0, 0.0, 0.0];
ϵ₀_2 = collect(-0.5:0.005:0.25);
type_2 = -1;
deg_2 = 1;

gap = ϵ₀_1 - ϵ₀_2;

μ = 0.0;
T = collect(50.:10:1000);

τ_form = Scattering.constant();             # relaxation time

# TENSORS COMPUTATION
σ = Array{Float64,2}(undef,length(T),length(gap));
S = Array{Float64,2}(undef,length(T),length(gap));
n = Array{Float64,2}(undef,length(T),length(gap));
# loop in which I move the cond band down and the val band up simultaneously
for i in 1:size(ϵ₀_1,1)
    band_1 = ParabBand(m_1,ϵ₀_1[i],type_1,deg_1);       # create the conduction band
    band_2 = ParabBand(m_2,ϵ₀_2[i],type_2,deg_2);       # create the valence band
    model = BandStructure(2,[band_1,band_2],μ);         # build the two-band structure
    σ[:,i] = electrical_conductivity(model,T,τ_form);
    S[:,i] = seebeck_coefficient(model,T,τ_form);
    n[:,i] = carrier_concentration(model,T,τ_form);
end

# PLOT
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

fig = plot(3, T, [σ,S*1e6,n], gap, titles=titles, xlabels=xlabels, ylabels=ylabels, zlabel=zlabel);

savefig(joinpath(savepath, "ex_two_isobands_convergence.png"), fig);
