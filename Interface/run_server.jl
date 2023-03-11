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


cd("../../Mstar2t")
using Pkg; Pkg.activate(".")

import HTTP
using Mstar2t: ComputingUnit


# create server with all the methods
const ROUTER = HTTP.Router()
HTTP.register!(ROUTER, "POST", "/api/clicalc", ComputingUnit.CLIcalc)
HTTP.register!(ROUTER, "POST", "/api/guicalc", ComputingUnit.GUIcalc)
HTTP.register!(ROUTER, "POST", "/api/guitaucalc", ComputingUnit.GUItaucalc)
HTTP.register!(ROUTER, "GET", "/api/check", ComputingUnit.check)

# run the server
HTTP.serve(ROUTER, "127.0.0.1", 1200, verbose=false)
