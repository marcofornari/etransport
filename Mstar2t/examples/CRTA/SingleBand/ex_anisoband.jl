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
import Mstar2t: Scattering
using LaTeXStrings

# BAND STRUCTURE DEFINITION
m = [0.1, 0.5, 3.0, 0.0, 0.0, 0.0];
ϵ₀ = 0.0;
type = -1;   # valence band
deg = 1;
band_1 = ParabBand(m,ϵ₀,type,deg);  # create the band

μ = collect(-.1:.01:.1);            # fermi level position
model = BandStructure(1,band_1,μ);  # build the band structure

T = collect(50:10:650);

τ_form = Scattering.constant();     # relaxation time

# TENSORS COMPUTATION
σ = electrical_conductivity(model,T,τ_form);
S = seebeck_coefficient(model,T,τ_form);
n = carrier_concentration(model,T,τ_form);

# PLOTS
titles = [L"$σ$ vs $T$, $τ = const$",
L"$S$ vs $T$, $τ = const$",
L"$n$ vs $T$, $τ = const$"];

xlabels = [L"$T$ $[K]$", 
L"$T$ $[K]$",
L"$T$ $[K]$"];

ylabels = [L"$\sigma$ $[(\Omega m)^{-1}]$",
L"$S$ $[\mu VK^{-1}]$",
L"n"];

zlabel = L"$μ$ $[eV]$";

fig = plot(3, T, [σ,S*1e6,n], μ, titles=titles, xlabels=xlabels, ylabels=ylabels, zlabel=zlabel);

savepath = joinpath("examples", "CRTA", "SingleBand");
savefig(joinpath(savepath, "ex_anisoband_val.png"), fig);
