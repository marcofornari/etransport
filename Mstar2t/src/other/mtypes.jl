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


@doc raw"""
    ParabBand(mstar::Vector{Float64},ϵ₀::Union{Float64,Vector{Float64}},type::Int64,deg::Int64)
 
Julia type to represent a band. Parameters:
mstar: effective mass tensor
ϵ₀: minimum/maximum of the band
type: conduction or valence band
deg: degeneracy
"""
struct ParabBand
    mstar::Vector{Float64}
    ϵ₀::Union{Float64,Vector{Float64}}
    type::Int64
    deg::Int64
end

@doc raw"""
    BandStructure(n::Int64, bands::Vector{ParabBand}, μ::Union{Float64,Vector{Float64}})
 
Julia type to represent a band full band structure. Parameters:
n: number of bands
bands: vector of bands
μ: Fermi level position
"""
struct BandStructure
    n::Int64
    bands::Vector{ParabBand}
    μ::Union{Float64,Vector{Float64}}
end

function BandStructure(n::Int64,bands::ParabBand,μ::Union{Float64,Vector{Float64}})
    return BandStructure(n,[bands],μ)
end

