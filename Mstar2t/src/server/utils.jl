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


# method to run simulations with CLI version of the interface of Mstar2t
function CLIcalc(req::HTTP.Request)
    result = runcalculation(JSON.parse(String(req.body)))
    if result == -10
        return HTTP.Response(210; body = JSON.json(result))
    elseif result == -20
        return HTTP.Response(210; body = JSON.json(result))
    elseif result == -30
        return HTTP.Response(210; body = JSON.json(result))
    end
    return HTTP.Response(200; body = JSON.json(result))
end


# method to run simulations with the GUI version of the interface of Mstar2t
function GUIcalc(req::HTTP.Request)
    result = runGUIcalculation(JSON.parse(String(req.body)))
    if result == -20
        return HTTP.Response(210; body = JSON.json(result)) # 210 code: relaxation time error
    elseif result == -30
        return HTTP.Response(210; body = JSON.json(result))
    elseif result == -40
        return HTTP.Response(210; body = JSON.json(result))
    end
    return HTTP.Response(200; body = JSON.json(result))
end


# method to compute tau 
function GUItaucalc(req::HTTP.Request)
    result = GUItaucalculation(JSON.parse(String(req.body)))
    if result == -20
        return HTTP.Response(210; body = JSON.json(result)) # 210 code: relaxation time error
    elseif result == -30
        return HTTP.Response(210; body = JSON.json(result))
    elseif result == -40
        return HTTP.Response(210; body = JSON.json(result))
    end
    return HTTP.Response(200; body = JSON.json(result))
end

# GUI method to check if server is up
function check(req::HTTP.Request)
    return HTTP.Response(200, JSON.json("server is up"))
end


# exception for Matthiessen's rule
struct MatthiessenError <: Exception end