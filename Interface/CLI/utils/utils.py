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



from PySide2 import QtCore, QtGui, QtWidgets


# show Mstar2t logo at launch
class LoadingScreen(QtWidgets.QMainWindow):
    def __init__(self, size):
        super().__init__()
        screen_width = size.width()
        screen_height = size.height()
        window_width = 280
        window_height = 160
        window_posx = (screen_width/2) - (window_width/2)
        window_posy = (screen_height/2) - (window_height/2)
        self.setWindowFlag(QtCore.Qt.FramelessWindowHint)
        self.setWindowTitle("")
        self.setGeometry(window_posx, window_posy, window_width, window_height)
        self.setWindowFlags(self.windowFlags() | QtCore.Qt.FramelessWindowHint)
        self.setAttribute(QtCore.Qt.WA_TranslucentBackground)
        self.setWindowIcon(QtGui.QIcon("icon.png"))
        self.label = QtWidgets.QLabel(self)
        self.label.setGeometry(QtCore.QRect(0, 0, 186, 106))
        self.label.setText("")
        self.label.setPixmap(QtGui.QPixmap("logo.png"))
        self.label.setScaledContents(True)
        self.label.setAlignment(QtCore.Qt.AlignCenter)
        self.label.setObjectName("label")
        QtCore.QTimer.singleShot(5000, self.close)
        self.show()
