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
import sys

import platform
import subprocess

from PySide2 import QtCore, QtWidgets

from utils.utils import LoadingScreen

if platform.system() == "Windows":
    import ctypes
    myappid = 'Mstar2t.bonal1l@cmich.edu' # arbitrary string
    ctypes.windll.shell32.SetCurrentProcessExplicitAppUserModelID(myappid)

# display
os.environ["QT_AUTO_SCREEN_SCALE_FACTOR"] = "1"
QtWidgets.QApplication.setAttribute(QtCore.Qt.AA_EnableHighDpiScaling, True)

# logo
app = QtWidgets.QApplication(sys.argv)
screen = app.primaryScreen()
loading = LoadingScreen(screen.size())
app.exec_()

# server
server = subprocess.Popen(['julia', '../run_server.jl'])
