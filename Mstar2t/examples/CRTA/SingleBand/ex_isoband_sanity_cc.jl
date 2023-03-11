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


# physics constants
const hbar = 6.62607e-34/(2*pi) # Reduced Planck's Constant
const qe = 1.602176565000000036063791900824543565343e-19 # Absolute charge of an electron (in Coulombs)
const me = 9.10938e-31 # Mass of an electron (in kg)
const kB = 1.38065e-23 # Boltzmann's Constant (In J/K)

# carrier concentration for a metal and insulator from statistical mechanics (grand canonical ensemble) 
# Ref: N. A. Mecholsky, L. Resca, I. L. Pegg, M. Fornari, Phys. Rev. B 2014, 89, 155131. 
function ngc_metal(band,ϵ,μ)
    return (-2*me*((band[1]*band[2]*band[3])^(1/3))*(ϵ-μ)*qe/((3^(2/3))*(hbar^2)*(pi^(4/3))))^(3/2)
end

function ngc_ins(band,ϵ,μ,t)
    β = /(1,(kB*t))
    return (me/(pi*β))^(3/2)*sqrt(band[1]*band[2]*band[3])*exp(β*(μ-ϵ)*qe)/(sqrt(2)*hbar^3)
end

# BAND STRUCTURE DEFINITION
m = [1.0, 1.0, 1.0, 0.0, 0.0, 0.0];
ϵ₀ = 0.0;
type = 1;   # conduction band
deg = 1;
band_1 = ParabBand(m,ϵ₀,type,deg);  # create the band

μ = collect(-1:0.5:3);            # fermi level position
model = BandStructure(1,band_1,μ)   # build the band structure
T = collect(50:10:650);
τ_form = Scattering.constant()      # relaxation time

# TENSORS COMPUTATION
σ = electrical_conductivity(model,T,τ_form);
S = seebeck_coefficient(model,T,τ_form);
n = carrier_concentration(model,T,τ_form);

# PLOT
titles = [L"$σ$ vs $T$, $τ = const$",
L"$S$ vs $T$, $τ = const$",
L"$n$ vs $T$, $τ = const$"]

xlabels = [L"$T$ $[K]$", 
L"$T$ $[K]$",
L"$T$ $[K]$"]

ylabels = [L"$\sigma$ $[(\Omega m)^{-1}]$",
L"$S$ $[\mu VK^{-1}]$",
L"$n$ $[m^{-3}]$"]

zlabel = L"$μ$ $[eV]$"

savepath = joinpath("examples", "CRTA", "SingleBand");

fig = plot(3, T, [σ,S*1e6,n], μ, titles=titles, xlabels=xlabels, ylabels=ylabels, zlabel=zlabel);

argexp = Array{Float64,2}(undef,length(μ),length(T));
n_ins = Array{Float64,2}(undef,length(μ),length(T));
for i in eachindex(T)
    for j in eachindex(μ)
        argexp[j,i] = 1/(kB*T[i])*(ϵ₀-μ[j])*qe
        n_ins[j,i] = ngc_ins(m,ϵ₀,μ[j],T[i])
    end 
end

CairoMakie.scatter!(fig[1,1], T, fill(qe*qe*1e-14/me*n[end],length(T)), markersize=3, marker='◀', color=:black)
CairoMakie.scatter!(fig[1,1], T, fill(qe*qe*1e-14/me*ngc_metal(m,ϵ₀,μ[end]),length(T)), markersize=3, color=:gray)
CairoMakie.scatter!(fig[1,3], T, fill(ngc_metal(m,ϵ₀,μ[end]),length(T)), markersize=3, color=:gray)
CairoMakie.lines!(fig[1,3], T, n_ins[1,:], linewidth=1, color=:red)
savefig(joinpath(savepath, "ex_isoband_sanity.png"), fig);