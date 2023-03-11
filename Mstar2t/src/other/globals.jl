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


# Physics constants
const hbar      = 6.62607e-34/(2*pi) # Reduced Planck's Constant
const qe        = 1.602176565000000036063791900824543565343e-19 # Absolute charge of an electron (in Coulombs)
const kB        = 1.38065e-23 # Boltzmann's Constant (In J/K)
const me        = 9.10938e-31 # Mass of an electron (in kg)
const mu0       = 1.602176e-19 # 1 ev in Joules (J/eV) 
const sl        = 299792458 # speed of light (m/s)
const ϵijk      = cat([0 0 0; 0 0 1; 0 -1 0], [0 0 -1; 0 0 0; 1 0 0], [0 1 0; -1 0 0; 0 0 0], dims=3)   # levi civta tensor in 3 dims

# Mstar2t parameters
const Nint      = 200 # Size of the Gaussian Quadrature grid

# precalculation of constants
const inv_kB    = 1/kB
const inv_mu0   = 1/qe 
const pref      = sqrt(me)/(2^(3/2)*pi^3*hbar^3)
const sspref    = qe^3/(me*sl)
const ccpref    = (2*me)^(3/2)/(3*pi^2*hbar^3)