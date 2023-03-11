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


module Transport


using LinearAlgebra
import SpecialFunctions
using FastGaussQuadrature
import HypergeometricFunctions
import PolyLog
using QuadGK
import Einsum
using SharedArrays

using CairoMakie: Figure

using ..Utils
using ..Scattering: ScModel, Matthiessen, compute_τ, matthiessen_rule

export  precalculation!,
        electrical_conductivity,
        seebeck_coefficient,
        carrier_concentration,
        thermal_conductivity,
        compute_elcond, compute_seebeck, 
        compute_carrierconc, compute_thermalcond,
        hall,
        lorenz_tensor,
        tint

include("../other/globals.jl")

include("_transport.jl")
include("integration.jl")

end