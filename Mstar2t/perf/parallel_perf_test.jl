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
using BenchmarkTools
using Mstar2t: Scattering

function test()
    # BAND STRUCTURE DEFINITION
    m_1 = [1.0, 1.0, 1.0, 0.0, 0.0, 0.0];
    ϵ₀_1 = 1.0;
    type_1 = 1;
    deg_1 = 1;
    band_1 = ParabBand(m_1,ϵ₀_1,type_1,deg_1);   # create the conduction band
    
    m_2 = [0.1, 1.0, 5.0, 0.0, 0.0, 0.0];
    ϵ₀_2 = 0.0;
    type_2 = -1;
    deg_2 = 1;
    band_2 = ParabBand(m_2,ϵ₀_2,type_2,deg_2);   # create the valence band

    μ = collect(-.1:0.01:.1);
    model = BandStructure(2,[band_1,band_2],μ);   # build the two-band structure

    T = collect(50.:10:750);

    τ_form = Scattering.constant();

    num_calc = length(μ)*length(T)
    println("Number of calculation for each tensor: ", num_calc)

    # TENSORS COMPUTATION
    @btime electrical_conductivity($model,$T,$τ_form);
    @btime seebeck_coefficient($model,$T,$τ_form);
    @btime carrier_concentration($model,$T,$τ_form);
end;

test();