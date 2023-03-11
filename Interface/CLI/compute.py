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


import os
import argparse
import subprocess

from colorama import Fore, Style

import requests
from datetime import datetime

from utils.reading_class import ReadInput

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--inputfile',
                    required=True,
                    help='path to input file')
parser.add_argument("--conductivity", "-e",
                    help="compute electrical conductivity",
                    action='store_true')
parser.add_argument("--seebeck", "-s",
                    help="compute Seebeck coefficient",
                    action='store_true')
parser.add_argument("--thermal", "-k",
                    help="compute thermal conductivity",
                    action='store_true')
parser.add_argument("--concentration", "-n",
                    help="compute carrier concentration",
                    action='store_true')
parser.add_argument("--tplot",
                    help="Temperature plot of the tensors",
                    action='store_true')
parser.add_argument("--muplot",
                    help="Fermi level plot of the results",
                    action='store_true')

args = parser.parse_args()
# get path of input file
data_path = args.inputfile

# 1. read the params from input file and create a python dict of parameters
params = ReadInput(data_path).read_params()

# 2. add the command line arguments to the python dict of parameters
dict_args = vars(args)
dict_args.pop("inputfile")
params["args"] = dict_args

# 3. Send a calculation request
url_c = 'http://127.0.0.1:1200/api/clicalc'
headers = {'Accept': 'application/json', 'Content-Type': 'application/json'}

try:
    r_calc = requests.post(url_c, json=params, headers=headers)
    r_calc.raise_for_status()
    # 4. Check response
    if r_calc.status_code == 200:
        return_value = r_calc.text
        print("Computation done.")
        print("Check results in " + return_value)
    elif r_calc.status_code == 210 and r_calc.text == "-10":
        print(f"{Fore.RED + Style.BRIGHT}ERROR:{Style.RESET_ALL} Export path not found.")
        print("Check input file [results path].")
    elif r_calc.status_code == 210 and r_calc.text == "-20":
        print(f"{Fore.RED + Style.BRIGHT}ERROR:{Style.RESET_ALL} Relaxation time functional form unknown.")
        print("Check input file.")
    elif r_calc.status_code == 210 and r_calc.text == "-30":
        print(f"{Fore.RED + Style.BRIGHT}ERROR:{Style.RESET_ALL} Domain error in the τ function calculation.")
        print("Possible reasons: negative arguments in logs, negative radicand in radical expressions.")
        print(f"{Fore.CYAN + Style.BRIGHT}Hint:{Style.RESET_ALL} shift the zero value of the chosen energy scale in the constraction of the band structure.")
    else:
        print(f"{Fore.RED + Style.BRIGHT}ERROR:{Style.RESET_ALL} Check server terminal.")
        print(r_calc)

except requests.exceptions.HTTPError as errh:
    print("Http error:", errh)
except requests.exceptions.ConnectionError as errc:
    print("Exception occurred (ConnectionError). Check server terminal.")
except requests.exceptions.Timeout as errt:
    print("Timeout error:", errt)
except requests.exceptions.RequestException as err:
    print("Exception occurred (RequestException). Check server terminal.", err)



