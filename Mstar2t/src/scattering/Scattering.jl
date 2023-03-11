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


module Scattering

using Parameters

using ..Utils: ParabBand, BandStructure, get_enlowestband

export  ScModel, Matthiessen,
        constant, T_fun, acoustic, impurity, matthiessen, 
        compute_τ, matthiessen_rule

        
abstract type ScModel end

@with_kw struct Acoustic <: ScModel
    τ::Function
    name::String = "acoustic"
end

@with_kw struct Impurity <: ScModel
    τ::Function
    name::String = "impurity"
end

@with_kw struct Constant <: ScModel
    τ::Function
    name::String = "constant"
end

@with_kw struct TFunc <: ScModel
    τ::Function
    name::String = "tfun"
end

@with_kw struct Matthiessen
    τ_models::Array{ScModel}
    γ::Real = -1
end


include("../other/globals.jl")
include("models.jl")

end