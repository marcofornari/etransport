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


module Plot

using CairoMakie
using PlotUtils
using LaTeXStrings
using DataFrames: DataFrame

using ..Utils: ParabBand, BandStructure, to_matrix
using ..Scattering: ScModel, Acoustic, Matthiessen, compute_τ, matthiessen_rule
using ..Transport: tint

export  plot,
        plot!,
        tplot_allcomp,
        plot_bandstructure,
        plot_τ,
        plot_bandcomp,
        savefig

include("_plot.jl")

end