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


import re


# this class takes care of reading and parsing the parameters of the 
# input file and storing them into python data structures. 
class ReadInput(object):

    def __init__(self, data_path):
        self.data_path = data_path
        self.num_bands = None
        # dict of input parameters
        self.data_input = {"# results fullpath": None,
                           "# export all data [true/false]": None,
                           "# number of bands": None,
                           "# Fermi level": None,
                           "# temperature": None,
                           "# bands masses and angles": None,
                           "# band type": None,
                           "# energy extrema": None,
                           "# degeneracy": None,
                           "# tau model [constant/acoustic/impurity]": None,
                           "# tau acoustic coefficients": None,
                           "# tau impurity coefficients": None }
        # dict of acoustic relaxation time params
        self.tauacoustic_coeffs = {"ϵ_min": 0.0,
                                   "A_sm": 1.0,
                                   "τm_max": 1.0,
                                   "T₀": 50.0,
                                   "μ_min": 2.0,
                                   "μ_max": 2.0}
        # dict of impurity relaxation time params
        self.tauimpurity_coeffs = {"ϵ_im": 1.0,
                                   "A_im": 1.0,
                                   "γ_im": 1.0}

    # function to add a new pair of key and value to data_input dict
    def add2dict(self, f):
        key = f.readline().rstrip('\n')
        value = f.readline().rstrip('\n')
        self.data_input[key] = value

    # function that reads the input file and fills the dictionary self.data_input above
    def read_params(self):
        with open(self.data_path, "r") as f:
            # output folder path -> self.data_input
            self.add2dict(f)
            # export_all_data variable -> self.data_input
            self.add2dict(f)
            # number of bands -> self.data_input
            key = f.readline().rstrip('\n')
            self.num_bands = int(f.readline().rstrip('\n'))
            self.data_input[key] = self.num_bands
            # temperature and fermi level -> self.data_input
            for i in range(2):
                self.add2dict(f)
            # bands masses, band type, energies and degeneracies (one value for each band) -> self.data_input
            for i in range(4):
                key = f.readline().rstrip('\n')
                value = list()
                for j in range(self.num_bands):
                    value.append(f.readline().rstrip('\n'))
                self.data_input[key] = value
            # tau model type -> self.data_input
            self.add2dict(f)
            # tau acoustic coefficients -> self.data_input
            key = f.readline().rstrip('\n')
            value = list()
            for j in range(len(self.tauacoustic_coeffs)):
                value.append(float(re.findall("\d*\.?\d+", f.readline().rstrip('\n'))[-1]))
            self.data_input[key] = value
            # tau impurity coefficients -> self.data_input
            key = f.readline().rstrip('\n')
            value = list()
            for j in range(len(self.tauimpurity_coeffs)):
                value.append(float(re.findall("\d*\.?\d+", f.readline().rstrip('\n'))[-1]))
            self.data_input[key] = value
        return self.data_input


