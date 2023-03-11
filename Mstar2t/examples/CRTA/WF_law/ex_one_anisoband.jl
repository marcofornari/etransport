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
using Colors


# BAND STRUCTURE DEFINITION
m = [0.1, 1.0, 10.0, 0.0, 0.0, 0.0];
ϵ₀ = 0.0;
type = 1;
deg = 1;
band_1 = ParabBand(m,ϵ₀,type,deg);   # create the band

μ = collect(-0.2:0.0005:0.2);
model = BandStructure(1,band_1,μ);   # build the band structure

T = collect(50.:50:1000);

τ_form = Scattering.constant();     # relaxation time

# TENSORS COMPUTATION
σ = electrical_conductivity(model,T,τ_form);
K = thermal_conductivity(model,T,τ_form);
L = lorenz_tensor(model,T,τ_form);

# PLOTS of the simulated tensors using Mstar2t
titles = [L"$σ$ vs $μ$, $τ = const$",
          L"$K$ vs $μ$, $τ = const$",
          L"$L$ vs $μ$, $τ = const$"];

xlabels = [L"$μ$ $[eV]$" for i in 1:3]; 

ylabels = [L"$\sigma$ $[(\Omega m)^{-1}]$",
L"$K$ $[W$ $(mK)^{-1}]$",
L"$L$ $[W \Omega K^{-2}]$"];

zlabel = L"$T$ $[K]$";

savepath = joinpath("examples", "CRTA", "WF_law");

fig = plot(3, μ, [Matrix(transpose(σ)),Matrix(transpose(K)),Matrix(transpose(L))], T, titles=titles, xlabels=xlabels, ylabels=ylabels, zlabel=zlabel);

# ANALYTIC solution
# The following lines can be used to plot the analytic solutions
const qe = 1.602176565000000036063791900824543565343e-19 # Absolute charge of an electron (in Coulombs)
const kB = 1.38065e-23 # Boltzmann's Constant (In J/K)

# Lorentz number - theoretical value for a metal
function L_teo(T,μ,ϵ)
    s = qe*(ϵ-μ)/(kB*T)
    return kB^2*pi^2/(3qe^2)*(1 -7(5pi^2-3)/(120s^2))
    #return kB^2*pi^2/(3qe^2)*(1 -7(5pi^2-3)/(120s^2) +(31pi^4-7pi^2)/(320s^4))
end;

# create colormap
function get_colors(row,min)
    color = RGB[]
    for i in 1:length(row)
        if row[i] > min
            push!(color,colorant"gray")
        else
            push!(color,colorant"white")
        end
    end
    return color
end;

# lorentz teoric for a metal
μ_start = findfirst(μ .>= 0.01);    # start the loop from μ=μ_start (we just want ~metallic regime)
T_stop = findfirst(T .>= 300);
L_met = Array{Float64,2}(undef,length(T[1:T_stop]),length(μ[μ_start:end]));
for i in 1:length(μ[μ_start:end])
    for j in 1:length(μ[1:T_stop])
        L_met[j,i] = L_teo(T[j],μ[i+μ_start-1],ϵ₀)
    end
end;
# the following lines are just to constraint the grey lines properly
metal_mask = ((qe./(kB*T[1:T_stop]))*transpose(ϵ₀.-μ[μ_start:end]) .< -1) .* (1 .- (L_met .< minimum(L)));
L_backgr = minimum(L)*(1 .- metal_mask);
L_met = (L_met .* metal_mask) .+ L_backgr ;
# plot the dash lines
for row in eachrow(L_met)
    color = get_colors(row, minimum(L))
    CairoMakie.lines!(fig[1,3], μ[μ_start:end], row, color=color, linestyle=:dash)
end;
CairoMakie.lines!(fig[1,3], μ, fill((kB*pi)^2/(3qe^2),length(μ)), color=:gray, linestyle=:dash, label="analytic (met)");
CairoMakie.lines!(fig[1,3], μ, fill(5*kB^2/(2qe^2),length(μ)), color=:gray, linestyle=:dot, label="analytic (ins)");
CairoMakie.axislegend(position=:lt);

# finally export the figure
savefig(joinpath(savepath, "ex_one_anisoband_mu.png"), fig);
