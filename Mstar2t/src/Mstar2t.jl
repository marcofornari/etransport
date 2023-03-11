# *************************************************************************** #
# *                                                                         * #
# *         Mstar2t - Central Michigan University University, 2023          * #
# *                                                                         * #
# *************************************************************************** #
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


module Mstar2t

# required packages
using CSV, DataFrames


# imports
include("other/Utils.jl")
using .Utils
export  ParabBand,
        BandStructure,
        savedata

include("scattering/Scattering.jl")
using .Scattering
export  ScModel, 
        constant, acoustic, impurity, T_fun,
        compute_τ, matthiessen

include("transport/Transport.jl")
using .Transport
export  electrical_conductivity,
        seebeck_coefficient,
        carrier_concentration,
        thermal_conductivity,
        lorenz_tensor,
        hall

include("plot/Plot.jl")
using .Plot
export  plot,
        plot!,
        tplot_allcomp,
        plot_bandstructure,
        plot_τ,
        savefig
        
include("server/ComputingUnit.jl")
using .ComputingUnit
export  CLIcalc, GUIcalc

end