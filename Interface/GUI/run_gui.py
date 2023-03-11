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
import time
import threading

import numpy as np
import pandas as pd
import itertools

import requests
from PySide2 import QtCore, QtGui, QtWidgets

import matplotlib
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
import matplotlib.cm     as cm
matplotlib.use('Qt5Agg')

from utils.utils import decimal_digits, LoadingScreen, ClickableLineEdit, QRoundProgressBar


############ DESIGN PARAMETERS ############
# matplotlib params
plt.rcParams['font.size'] = 9
plt.rcParams['axes.titlesize'] = 14
plt.rcParams['axes.linewidth'] = 0.5
plt.rcParams['axes.labelsize'] = 9
plt.rcParams['axes.titlepad'] = 11
plt.rcParams['axes.spines.bottom'] = True
plt.rcParams['axes.spines.left'] = True
plt.rcParams['axes.spines.top'] = False
plt.rcParams['axes.spines.right'] = False
plt.rcParams['xtick.major.size'] = 2
plt.rcParams['xtick.labelsize'] = 7
plt.rcParams['ytick.major.size'] = 2
plt.rcParams['ytick.labelsize'] = 7
plt.rcParams['legend.fontsize'] = 7
plt.rcParams["legend.loc"] = "upper left"
plt.rcParams["lines.linewidth"] = 1
plt.rcParams['lines.markersize'] = 5

# general design params
gui_color = (31,161,135)
Qt_gui_color = QtGui.QColor(gui_color[0],gui_color[1],gui_color[2])
header_color = "rgb({}, {}, {})".format(gui_color[0],gui_color[1],gui_color[2])
mySetStyleSheet = """QHeaderView::section{
border-top:0px solid #D8D8D8;
border-left:0px solid #D8D8D8;
border-right:1px solid #D8D8D8;
border-bottom:1px solid #D8D8D8;
background-color: %s;
color:white;
padding:4px;}
QTableCornerButton::section{
border-top:0px solid #D8D8D8;
border-left:0px solid #D8D8D8;
border-right:1px solid #D8D8D8;
border-bottom: 1px solid #D8D8D8;
background-color:white;}""" % (header_color)

###########################################


############# DATA STRUCTURES #############
# class to store input params 
class Data(object):
    def __init__(self):
        self.export_all_data = False # TODO: implement button
        self.num_cbands = 0
        self.num_vbands = 0
        self.num_bands = 0
        self.cTensors = list()
        self.mband = list()
        self.ebandmins = list()
        self.degeneracies = list()
        self.T_str = "[300:1000:100]"
        self.mu_str = "0.5"
        self.T = 0
        self.mu = 0
        self.tau_model_type = "constant"
        self.tau_acoustic_coeffs = list()
        self.tau_impurity_coeffs = list()
        self.tau_matthiessen_models = "000"
        self.tau_matthiessen_gamma = "0.0"

# class to store experimental data
class ExpDataDB(object):
    def __init__(self):
        self.df = None
        self.verHeader = ['temperature', 'conductivity', 'seebeck', 'thermal', 'concentration']

    def __getitem__(self, key):
        return self.df[self.df.iloc[:, 0] == key].iloc[0, 1:].to_numpy()

    def set_exp_data(self, df):
        self.df = df
        rows = set(self.df.iloc[:, 0])
        for r in rows:
            if r not in self.verHeader:
                raise ValueError

    def clear(self):
        self.df = None


# class to store output data
class ResultAllCompData(object):
    def __init__(self):
        self.num_t = None
        self.num_mu = None
        self.conductivity = None
        self.seebeck = None
        self.thermal = None
        self.concentration = None

    def setTmu(self, num_mu, num_t):
        self.num_mu = num_mu
        self.num_t = num_t
        self.preallocate_tensors()

    def preallocate_tensors(self):
        self.conductivity = np.empty((6, self.num_mu, self.num_t))
        self.seebeck = np.empty((6, self.num_mu, self.num_t))
        self.thermal = np.empty((6, self.num_mu, self.num_t))
        self.concentration = np.empty((1, self.num_mu, self.num_t))

    def setCond(self, tensor):
        self.conductivity = tensor

    def setSeebeck(self, tensor):
        self.seebeck = tensor

    def setThermal(self, tensor):
        self.thermal = tensor

    def setConc(self, tensor):
        self.concentration = tensor

    def clear(self):
        self.conductivity = None
        self.seebeck = None
        self.thermal = None
        self.concentration = None

    def isallocated(self):
        return True if ((self.num_mu is not None) or (self.num_t is not None)) else False


# class to store output data (trace of each matrix)
class ResultTraceData(object):
    def __init__(self):
        self.data = []
        self.label = []

    def clear(self):
        self.data = []
        self.label = []

###########################################

# pyqt main window
class MainWindow(QtWidgets.QMainWindow):
    def __init__(self, server=None):
        super(MainWindow, self).__init__()
        self.server = server

    def closeEvent(self, event):
        if self.server is not None:
            self.server.terminate()
        event.accept()

# input window
class UiInputWindow(object):
    def __init__(self):
        self.is_first_run_thread_active = True
        self.is_first_run_completed = False
        self.data = Data()
        self.exp_data = ExpDataDB()
        self.out_trace_data = ResultTraceData()
        self.out_all_data = ResultAllCompData()

    def setupUi(self, InputWindow):
        self.InputWindow = InputWindow
        self.InputWindow.setObjectName("InputWindow")
        self.InputWindow.resize(730, 642)
        self.InputWindow.setAutoFillBackground(False)
        self.appIcon = QtGui.QIcon("./images/icon.png")
        self.InputWindow.setWindowIcon(self.appIcon)
        self.OutputWindow = QtWidgets.QMainWindow()
        self.ui_out = UiOutputWindow(self)
        self.ui_out.setupUi(self.OutputWindow)
        self.OutputWindow.show()
        self.centralwidget = QtWidgets.QWidget(self.InputWindow)
        self.centralwidget.setObjectName("centralwidget")
        self.gridFrame = QtWidgets.QFrame(self.centralwidget)
        self.gridFrame.setGeometry(QtCore.QRect(280, 5, 211, 41))
        self.gridFrame.setObjectName("gridFrame")
        self.gridLayout = QtWidgets.QGridLayout(self.gridFrame)
        self.gridLayout.setContentsMargins(-1, 1, -1, -1)
        self.gridLayout.setObjectName("gridLayout")
        self.InputWindow.setCentralWidget(self.centralwidget)

        # menubar
        self.menubar = QtWidgets.QMenuBar(self.InputWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 819, 20))
        self.menubar.setObjectName("menubar")
        self.menuFile = QtWidgets.QMenu(self.menubar)
        self.menuFile.setObjectName("menuFile")
        self.menuHelp = QtWidgets.QMenu(self.menubar)
        self.menuHelp.setObjectName("menuHelp")
        self.InputWindow.setMenuBar(self.menubar)
        # status bar
        self.statusbar = QtWidgets.QStatusBar(self.InputWindow)
        self.statusbar.setObjectName("statusbar")
        self.InputWindow.setStatusBar(self.statusbar)

        # bars actions
        self.actionImportData = QtWidgets.QAction(self.InputWindow)
        self.actionImportData.setObjectName("actionImportData")
        self.actionImportData.triggered.connect(self.import_exp_data)
        self.actionExit = QtWidgets.QAction(self.InputWindow)
        self.actionExit.setObjectName("actionExit")
        self.actionExit.triggered.connect(self.exit)
        self.actionAbout = QtWidgets.QAction(self.InputWindow)
        self.actionAbout.setObjectName("actionAbout")
        self.menubar.addAction(self.menuFile.menuAction())
        self.menubar.addAction(self.menuHelp.menuAction())
        self.menuFile.addAction(self.actionImportData)
        self.menuFile.addAction(self.actionExit)
        self.menuHelp.addAction(self.actionAbout)

        # header
        ## title
        font_title = QtGui.QFont()
        font_title.setPointSize(12)
        font_title.setBold(True)
        self.inputlabel = QtWidgets.QLabel(self.centralwidget)
        self.inputlabel.setGeometry(QtCore.QRect(269, 20, 193, 40))
        self.inputlabel.setFont(font_title)
        self.inputlabel.setTextFormat(QtCore.Qt.PlainText)
        self.inputlabel.setAlignment(QtCore.Qt.AlignCenter)
        self.inputlabel.setIndent(0)
        self.inputlabel.setObjectName("inputlabel")
        ## server status
        ### server box
        self.serverBox = QtWidgets.QGroupBox(self.centralwidget)
        self.serverBox.setGeometry(QtCore.QRect(537, 25, 150, 60))
        self.serverBox.setStyleSheet("QGroupBox#serverBox {background-color: none; border:0px;}")
        self.serverBox.setTitle("")
        self.serverBox.setObjectName("serverBox")
        ### server label
        font = QtGui.QFont()
        font.setPointSize(9)
        font.setBold(True)
        self.serverLabel = QtWidgets.QLabel(self.serverBox)
        self.serverLabel.setGeometry(QtCore.QRect(0, 5, 90, 20))
        self.serverLabel.setFont(font)
        self.serverLabel.setObjectName("serverLabel")
        ### down signal
        self.redSignal = QtWidgets.QLabel(self.serverBox)
        self.redSignal.setGeometry(QtCore.QRect(90, 5, 32, 20))
        self.redSignal.setScaledContents(True)
        self.redSignal.setPixmap(QtGui.QPixmap('./images/down.png'))
        self.redSignal.setObjectName('redSignal')
        ### up signal
        self.greenSignal = QtWidgets.QLabel(self.serverBox)
        self.greenSignal.setGeometry(QtCore.QRect(110, 5, 32, 20))
        self.greenSignal.setScaledContents(True)
        self.greenSignal.setPixmap(QtGui.QPixmap('./images/up.png'))
        self.greenSignal.setVisible(False)
        self.greenSignal.setObjectName('greenSignal')

        # INPUT section
        ## input box
        self.inputBox = QtWidgets.QGroupBox(self.centralwidget)
        self.inputBox.setGeometry(QtCore.QRect(55, 60, 620, 435))
        self.inputBox.setStyleSheet("QGroupBox{background-color: white}")
        self.inputBox.setTitle("")
        self.inputBox.setObjectName("inputGroupBox")

        ## conduction bands
        ### title
        self.CondText = QtWidgets.QLabel(self.inputBox)
        self.CondText.setGeometry(QtCore.QRect(9, 10, 235, 25))
        self.CondText.setStyleSheet("Background-color:white;border: 0.8px solid black;")
        font = QtGui.QFont()
        font.setPointSize(9)
        font.setBold(True)
        self.CondText.setFont(font)
        self.CondText.setObjectName("CondText")
        ### font
        ctfont = QtGui.QFont()
        ctfont.setBold(False)
        ctfont.setPointSize(8)
        ### spin box -> number of conduction bands
        self.CondSpin = QtWidgets.QSpinBox(self.inputBox)
        self.CondSpin.setGeometry(QtCore.QRect(243, 10, 62, 25))
        self.CondSpin.setObjectName("CondSpin")
        self.CondSpin.setValue(0)
        self.CondSpin.setFont(ctfont)
        self.CondSpin.valueChanged.connect(self.add_cond_band)
        ### conduction bands params table
        self.CondTable = QtWidgets.QTableWidget(self.inputBox)
        self.CondTable.setGeometry(QtCore.QRect(8, 33, 297, 120))
        self.CondTable.setMinimumSize(QtCore.QSize(201, 71))
        self.CondTable.setFont(ctfont)
        self.CondTable.setContextMenuPolicy(QtCore.Qt.ActionsContextMenu)
        self.CondTable.setAutoFillBackground(False)
        self.CondTable.setFrameShape(QtWidgets.QFrame.Box)
        self.CondTable.setFrameShadow(QtWidgets.QFrame.Raised)
        self.CondTable.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        self.CondTable.setLineWidth(1)
        self.CondTable.setMidLineWidth(1)
        self.CondTable.setTabKeyNavigation(True)
        self.CondTable.setDragDropOverwriteMode(False)
        self.CondTable.setAlternatingRowColors(True)
        self.CondTable.setVerticalScrollMode(QtWidgets.QAbstractItemView.ScrollPerPixel)
        self.CondTable.setSelectionMode(QtWidgets.QAbstractItemView.NoSelection)
        self.CondTable.setIconSize(QtCore.QSize(1, 1))
        self.CondTable.setShowGrid(True)
        self.CondTable.setGridStyle(QtCore.Qt.SolidLine)
        self.CondTable.setRowCount(0)
        self.CondTable.setColumnCount(4)
        item = QtWidgets.QTableWidgetItem()
        item.setTextAlignment(QtCore.Qt.AlignCenter)
        brush = QtGui.QBrush(QtGui.QColor(0, 0, 255))
        brush.setStyle(QtCore.Qt.NoBrush)
        item.setForeground(brush)
        self.CondTable.setVerticalHeaderItem(0, item)
        item = QtWidgets.QTableWidgetItem()
        self.CondTable.setVerticalHeaderItem(1, item)
        font_hh = QtGui.QFont()
        font_hh.setPointSize(8)
        font_hh.setBold(True)
        item = QtWidgets.QTableWidgetItem()
        item.setTextAlignment(QtCore.Qt.AlignCenter)
        item.setFont(font_hh)
        self.CondTable.setHorizontalHeaderItem(0, item)
        item = QtWidgets.QTableWidgetItem()
        item.setTextAlignment(QtCore.Qt.AlignCenter)
        item.setFont(font_hh)
        self.CondTable.setHorizontalHeaderItem(1, item)
        item = QtWidgets.QTableWidgetItem()
        item.setTextAlignment(QtCore.Qt.AlignCenter)
        item.setFont(font_hh)
        self.CondTable.setHorizontalHeaderItem(2, item)
        item = QtWidgets.QTableWidgetItem()
        item.setTextAlignment(QtCore.Qt.AlignCenter)
        item.setFont(font_hh)
        self.CondTable.setHorizontalHeaderItem(3, item)
        item = QtWidgets.QTableWidgetItem()
        item.setTextAlignment(QtCore.Qt.AlignCenter)
        self.CondTable.setItem(0, 0, item)
        item = QtWidgets.QTableWidgetItem()
        item.setTextAlignment(QtCore.Qt.AlignCenter)
        self.CondTable.setItem(0, 1, item)
        item = QtWidgets.QTableWidgetItem()
        item.setTextAlignment(QtCore.Qt.AlignCenter)
        self.CondTable.setItem(0, 2, item)
        item = QtWidgets.QTableWidgetItem()
        item.setTextAlignment(QtCore.Qt.AlignCenter)
        self.CondTable.setItem(0, 3, item)
        self.CondTable.horizontalHeader().setVisible(True)
        self.CondTable.horizontalHeader().setCascadingSectionResizes(True)
        self.CondTable.horizontalHeader().setDefaultSectionSize(73)
        self.CondTable.horizontalHeader().setHighlightSections(True)
        self.CondTable.horizontalHeader().setMinimumSectionSize(49)
        self.CondTable.horizontalHeader().setSortIndicatorShown(True)
        self.CondTable.horizontalHeader().setStretchLastSection(False)
        self.CondTable.horizontalHeader().setFixedHeight(28)
        self.CondTable.verticalHeader().setVisible(False)
        self.CondTable.verticalHeader().setHighlightSections(True)
        self.CondTable.verticalHeader().setSortIndicatorShown(False)
        self.CondTable.verticalHeader().setStretchLastSection(False)
        self.CondTable.setObjectName("CondTable")

        ## valence bands
        ### title
        self.ValText = QtWidgets.QLabel(self.inputBox)
        self.ValText.setGeometry(QtCore.QRect(314, 10, 245, 25))
        self.ValText.setStyleSheet("Background-color:white;border: 0.8px solid black;")
        font = QtGui.QFont()
        font.setPointSize(9)
        font.setBold(True)
        self.ValText.setFont(font)
        self.ValText.setObjectName("ValText")
        ### font
        vtfont = QtGui.QFont()
        vtfont.setPointSize(8)
        vtfont.setBold(False)
        ### spin box -> number of valence bands
        self.ValSpin = QtWidgets.QSpinBox(self.inputBox)
        self.ValSpin.setGeometry(QtCore.QRect(550, 10, 60, 25))
        self.ValSpin.setObjectName("ValSpin")
        self.ValSpin.setValue(0)
        self.ValSpin.setFont(vtfont)
        self.ValSpin.valueChanged.connect(self.add_val_band)
        ### valence bands params table
        self.ValTable = QtWidgets.QTableWidget(self.inputBox)
        self.ValTable.setGeometry(QtCore.QRect(313, 33, 297, 120))
        self.ValTable.setMinimumSize(QtCore.QSize(201, 72))
        self.ValTable.setFont(vtfont)
        self.ValTable.setMouseTracking(False)
        self.ValTable.setContextMenuPolicy(QtCore.Qt.ActionsContextMenu)
        self.ValTable.setAutoFillBackground(False)
        self.ValTable.setFrameShape(QtWidgets.QFrame.Box)
        self.ValTable.setFrameShadow(QtWidgets.QFrame.Raised)
        self.ValTable.setLineWidth(1)
        self.ValTable.setMidLineWidth(1)
        self.ValTable.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        self.ValTable.setTabKeyNavigation(True)
        self.ValTable.setDragDropOverwriteMode(False)
        self.ValTable.setAlternatingRowColors(True)
        self.ValTable.setVerticalScrollMode(QtWidgets.QAbstractItemView.ScrollPerPixel)
        self.ValTable.setSelectionMode(QtWidgets.QAbstractItemView.NoSelection)
        self.ValTable.setIconSize(QtCore.QSize(1, 1))
        self.ValTable.setShowGrid(True)
        self.ValTable.setGridStyle(QtCore.Qt.SolidLine)
        self.ValTable.setRowCount(0)
        self.ValTable.setColumnCount(4)
        self.ValTable.horizontalHeader().setFont(font)
        font = self.ValTable.horizontalHeader().font()
        item = QtWidgets.QTableWidgetItem()
        item.setTextAlignment(QtCore.Qt.AlignCenter)
        item.setFont(font_hh)
        self.ValTable.setVerticalHeaderItem(0, item)
        item = QtWidgets.QTableWidgetItem()
        item.setTextAlignment(QtCore.Qt.AlignCenter)
        item.setFont(font_hh)
        self.ValTable.setHorizontalHeaderItem(0, item)
        item = QtWidgets.QTableWidgetItem()
        item.setFont(font_hh)
        item.setTextAlignment(QtCore.Qt.AlignCenter)
        self.ValTable.setHorizontalHeaderItem(1, item)
        item = QtWidgets.QTableWidgetItem()
        item.setTextAlignment(QtCore.Qt.AlignCenter)
        item.setFont(font_hh)
        self.ValTable.setHorizontalHeaderItem(2, item)
        item = QtWidgets.QTableWidgetItem()
        item.setTextAlignment(QtCore.Qt.AlignCenter)
        item.setFont(font_hh)
        self.ValTable.setHorizontalHeaderItem(3, item)
        item = QtWidgets.QTableWidgetItem()
        item.setTextAlignment(QtCore.Qt.AlignCenter)
        self.ValTable.setItem(0, 0, item)
        item = QtWidgets.QTableWidgetItem()
        item.setTextAlignment(QtCore.Qt.AlignCenter)
        self.ValTable.setItem(0, 1, item)
        item = QtWidgets.QTableWidgetItem()
        item.setTextAlignment(QtCore.Qt.AlignCenter)
        self.ValTable.setItem(0, 2, item)
        item = QtWidgets.QTableWidgetItem()
        item.setTextAlignment(QtCore.Qt.AlignCenter)
        self.ValTable.setItem(0, 3, item)
        self.ValTable.horizontalHeader().setVisible(True)
        self.ValTable.horizontalHeader().setCascadingSectionResizes(True)
        self.ValTable.horizontalHeader().setDefaultSectionSize(73)
        self.ValTable.horizontalHeader().setHighlightSections(True)
        self.ValTable.horizontalHeader().setMinimumSectionSize(49)
        self.ValTable.horizontalHeader().setSortIndicatorShown(True)
        self.ValTable.horizontalHeader().setStretchLastSection(False)
        self.ValTable.horizontalHeader().setFixedHeight(25)
        self.ValTable.verticalHeader().setVisible(False)
        self.ValTable.verticalHeader().setHighlightSections(True)
        self.ValTable.verticalHeader().setSortIndicatorShown(False)
        self.ValTable.verticalHeader().setStretchLastSection(False)
        self.ValTable.setSelectionMode(QtWidgets.QAbstractItemView.NoSelection)
        self.ValTable.setObjectName("ValTable")

        ## Thermodynamics params
        self.tmuBox = QtWidgets.QGroupBox(self.inputBox)
        self.tmuBox.setGeometry(QtCore.QRect(10, 205, 200, 54))
        self.tmuBox.setStyleSheet("Background-color:white")
        self.tmuBox.setTitle("")
        self.tmuBox.setObjectName("tmuBox")
        ### temperature title
        self.tLabel = QtWidgets.QLabel(self.tmuBox)
        self.tLabel.setGeometry(QtCore.QRect(0, 0, 100, 28))
        self.tLabel.setStyleSheet("Background-color: {}; border-top: 2px solid darkGray; border-top: 2px solid darkGray; border-bottom: 2px solid darkGray; border-left: 2px solid darkGray; border-right: 1px solid lightGray; color:white".format(header_color))
        self.tLabel.setFont(font_hh)
        self.tLabel.setTextFormat(QtCore.Qt.PlainText)
        self.tLabel.setAlignment(QtCore.Qt.AlignCenter)
        self.tLabel.setIndent(0)
        self.tLabel.setObjectName("tLabel")
        ### temperature input 
        font = QtGui.QFont()
        font.setPointSize(9)
        self.tInput = ClickableLineEdit(self.tmuBox)
        self.tInput.setGeometry(QtCore.QRect(0, 25, 100, 28))
        self.tInput.setStyleSheet("QLineEdit {Background-color:white; border-top: 2px solid darkGray; border-top: 2px solid darkGray; border-bottom: 2px solid darkGray; border-left: 2px solid darkGray; border-right: 1px solid lightGray; color: gray}")        
        self.tInput.clicked.connect(self.insert_new_temperature)
        self.tInput.editingFinished.connect(self.update_temp)
        self.tInput.setFont(font)
        self.tInput.setAlignment(QtCore.Qt.AlignCenter)
        self.tInput.setObjectName("tInput")
        ### Fermi level title
        self.muLabel = QtWidgets.QLabel(self.tmuBox)
        self.muLabel.setGeometry(QtCore.QRect(99, 0, 100, 28))
        self.muLabel.setStyleSheet("Background-color: {}; border-top: 2px solid darkGray; border-top: 2px solid darkGray; border-bottom: 2px solid darkGray; border-right: 2px solid darkGray; border-left: 1px solid lightGray; color:white".format(header_color))
        self.muLabel.setFont(font_hh)
        self.muLabel.setTextFormat(QtCore.Qt.PlainText)
        self.muLabel.setAlignment(QtCore.Qt.AlignCenter)
        self.muLabel.setIndent(0)
        self.muLabel.setObjectName("muLabel")
        ### Fermi level input
        self.muInput = ClickableLineEdit(self.tmuBox)
        self.muInput.setGeometry(QtCore.QRect(99, 25, 100, 28))
        self.muInput.setStyleSheet("QLineEdit {Background-color:white; border-top: 2px solid darkGray; border-top: 2px solid darkGray; border-bottom: 2px solid darkGray; border-right: 2px solid darkGray; border-left: 1px solid lightGray; color: gray}")
        self.muInput.clicked.connect(self.insert_new_mu)
        self.muInput.editingFinished.connect(self.update_mu)
        self.muInput.setFont(font)
        self.muInput.setAlignment(QtCore.Qt.AlignCenter)
        self.muInput.setObjectName("muInput")

        ## Relexation time box
        self.tauBox = QtWidgets.QGroupBox(self.inputBox)
        self.tauBox.setGeometry(QtCore.QRect(233, 160, 377, 265))
        self.tauBox.setStyleSheet("background-color: white; border:1px solid darkGray;")
        self.tauBox.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.tauBox.setObjectName("tauBox")
        ### relexation time title
        self.tauTitle = QtWidgets.QLabel(self.tauBox)
        self.tauTitle.setGeometry(QtCore.QRect(143, 5, 80, 20))
        font = QtGui.QFont()
        font.setPointSize(9)
        font.setBold(True)
        self.tauTitle.setStyleSheet("background-color: white; border:0px;")
        self.tauTitle.setFont(font)
        self.tauTitle.setObjectName("tauTitle")
        ### relexation time combobox [constant/acoustic/impurity/matthiessen]
        self.comboBox = QtWidgets.QComboBox(self.tauBox)
        self.comboBox.setGeometry(QtCore.QRect(204, 5, 115, 20))
        font = QtGui.QFont()
        font.setPointSize(9)
        self.comboBox.setFont(font)
        self.comboBox.addItem("")
        self.comboBox.addItem("")
        self.comboBox.addItem("")
        self.comboBox.addItem("")
        # self.comboBox.currentIndexChanged.connect(self.computetau)
        self.comboBox.setObjectName("comboBox")
        ### relexation time plot button
        self.TauplotButton = QtWidgets.QPushButton(self.tauBox)
        self.TauplotButton.setGeometry(QtCore.QRect(322, 5, 50, 20))
        font = QtGui.QFont()
        font.setPointSize(8)
        self.TauplotButton.setFont(font)
        self.TauplotButton.setLayoutDirection(QtCore.Qt.LeftToRight)
        self.TauplotButton.setAutoFillBackground(True)
        self.TauplotButton.setStyleSheet("QPushButton {background-color: lightgray;} QPushButton:hover{background-color: %s;}" % (header_color))
        self.TauplotButton.blockSignals(True)
        self.TauplotButton.clicked.connect(self.computetau)
        self.TauplotButton.setObjectName("TauplotButton")
        ### acoustic relexation time coeff
        #### title
        self.acTitle = QtWidgets.QLabel(self.tauBox)
        self.acTitle.setGeometry(QtCore.QRect(295, 28, 75, 20))
        font = QtGui.QFont()
        font.setPointSize(7)
        self.acTitle.setStyleSheet("background-color: white; border:0px;")
        self.acTitle.setFont(font)
        self.acTitle.setObjectName("acTitle")
        #### table
        font_ac = QtGui.QFont()
        font_ac.setPointSize(8)
        self.TauAcousticCoeffTable = QtWidgets.QTableWidget(self.tauBox)
        self.TauAcousticCoeffTable.setStyleSheet("background-color: white; border:1px solid darkGray;")
        self.TauAcousticCoeffTable.setGeometry(QtCore.QRect(5, 45, 367, 52))
        self.TauAcousticCoeffTable.setMinimumSize(QtCore.QSize(0, 0))
        self.TauAcousticCoeffTable.setFont(font_ac)
        self.TauAcousticCoeffTable.setFrameShape(QtWidgets.QFrame.Box)
        self.TauAcousticCoeffTable.setFrameShadow(QtWidgets.QFrame.Raised)
        self.TauAcousticCoeffTable.setLineWidth(1)
        self.TauAcousticCoeffTable.setMidLineWidth(1)
        self.TauAcousticCoeffTable.setAutoFillBackground(False)
        self.TauAcousticCoeffTable.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        self.TauAcousticCoeffTable.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        self.TauAcousticCoeffTable.setSelectionMode(QtWidgets.QAbstractItemView.NoSelection)
        self.TauAcousticCoeffTable.setShowGrid(True)
        self.TauAcousticCoeffTable.setGridStyle(QtCore.Qt.SolidLine)
        self.TauAcousticCoeffTable.setCornerButtonEnabled(True)
        self.TauAcousticCoeffTable.setColumnCount(6)
        self.TauAcousticCoeffTable.setRowCount(1)
        item = QtWidgets.QTableWidgetItem()
        self.TauAcousticCoeffTable.setVerticalHeaderItem(0, item)
        font_tau = QtGui.QFont()
        font_tau.setPointSize(8)
        font_tau.setBold(True)
        item = QtWidgets.QTableWidgetItem()
        item.setTextAlignment(QtCore.Qt.AlignCenter)
        item.setFont(font_tau)
        self.TauAcousticCoeffTable.setHorizontalHeaderItem(0, item)
        item = QtWidgets.QTableWidgetItem()
        item.setTextAlignment(QtCore.Qt.AlignCenter)
        item.setFont(font_tau)
        self.TauAcousticCoeffTable.setHorizontalHeaderItem(1, item)
        item = QtWidgets.QTableWidgetItem()
        item.setTextAlignment(QtCore.Qt.AlignCenter)
        item.setFont(font_tau)
        self.TauAcousticCoeffTable.setHorizontalHeaderItem(2, item)
        item = QtWidgets.QTableWidgetItem()
        item.setTextAlignment(QtCore.Qt.AlignCenter)
        item.setFont(font_tau)
        self.TauAcousticCoeffTable.setHorizontalHeaderItem(3, item)
        item = QtWidgets.QTableWidgetItem()
        item.setTextAlignment(QtCore.Qt.AlignCenter)
        item.setFont(font_tau)
        self.TauAcousticCoeffTable.setHorizontalHeaderItem(4, item)
        item = QtWidgets.QTableWidgetItem()
        item.setTextAlignment(QtCore.Qt.AlignCenter)
        item.setFont(font_tau)
        self.TauAcousticCoeffTable.setHorizontalHeaderItem(5, item)
        item = QtWidgets.QTableWidgetItem()
        item.setTextAlignment(QtCore.Qt.AlignCenter)
        self.TauAcousticCoeffTable.setItem(0, 0, item)
        item = QtWidgets.QTableWidgetItem()
        item.setTextAlignment(QtCore.Qt.AlignCenter)
        self.TauAcousticCoeffTable.setItem(0, 1, item)
        item = QtWidgets.QTableWidgetItem()
        item.setTextAlignment(QtCore.Qt.AlignCenter)
        self.TauAcousticCoeffTable.setItem(0, 2, item)
        item = QtWidgets.QTableWidgetItem()
        item.setTextAlignment(QtCore.Qt.AlignCenter)
        self.TauAcousticCoeffTable.setItem(0, 3, item)
        item = QtWidgets.QTableWidgetItem()
        item.setTextAlignment(QtCore.Qt.AlignCenter)
        self.TauAcousticCoeffTable.setItem(0, 4, item)
        item = QtWidgets.QTableWidgetItem()
        item.setTextAlignment(QtCore.Qt.AlignCenter)
        self.TauAcousticCoeffTable.setItem(0, 5, item)
        self.TauAcousticCoeffTable.horizontalHeader().setVisible(True)
        self.TauAcousticCoeffTable.horizontalHeader().setCascadingSectionResizes(True)
        self.TauAcousticCoeffTable.horizontalHeader().setDefaultSectionSize(15)
        self.TauAcousticCoeffTable.horizontalHeader().setHighlightSections(True)
        self.TauAcousticCoeffTable.horizontalHeader().setMinimumSectionSize(15)
        self.TauAcousticCoeffTable.verticalHeader().setVisible(False)
        self.TauAcousticCoeffTable.horizontalHeader().setFixedHeight(25)
        self.TauAcousticCoeffTable.verticalHeader().setSortIndicatorShown(False)
        self.TauAcousticCoeffTable.verticalHeader().setStretchLastSection(False)
        self.TauAcousticCoeffTable.horizontalHeader().setSectionResizeMode(QtWidgets.QHeaderView.Stretch)
        self.TauAcousticCoeffTable.verticalHeader().setSectionResizeMode(QtWidgets.QHeaderView.Stretch)
        self.TauAcousticCoeffTable.setSizeAdjustPolicy(QtWidgets.QAbstractScrollArea.AdjustToContents)
        self.TauAcousticCoeffTable.itemChanged.connect(self.computetau)
        self.TauAcousticCoeffTable.setObjectName("TauAcousticCoeffTable")
        #### acoustic relaxation time e_min param box
        self.acousBox = QtWidgets.QComboBox(self.tauBox)
        self.acousBox.setGeometry(QtCore.QRect(5, 98, 120, 20))
        font = QtGui.QFont()
        font.setPointSize(9)
        self.acousBox.setFont(font)
        self.acousBox.addItem("")
        self.acousBox.setObjectName("acousBox")
        #### update acousbox button
        self.updateAcousBoxButton = QtWidgets.QPushButton(self.tauBox)
        self.updateAcousBoxButton.setStyleSheet("QPushButton#updateAcousBoxButton {border-image: url(./images/update.png); background-repeat: no-repeat;} QPushButton#updateAcousBoxButton:hover{border-image: url(./images/updatehover.png); background-repeat: no-repeat;} ")
        self.updateAcousBoxButton.setGeometry(128, 98, 20, 20)
        self.updateAcousBoxButton.setObjectName("updateAcousBoxButton")
        self.updateAcousBoxButton.clicked.connect(self.update_ebandmins)
        ### impurity relexation time coeff
        #### title
        self.imTitle = QtWidgets.QLabel(self.tauBox)
        self.imTitle.setGeometry(QtCore.QRect(294, 105, 75, 15))
        font = QtGui.QFont()
        font.setPointSize(7)
        self.imTitle.setStyleSheet("background-color: white; border:0px;")
        self.imTitle.setFont(font)
        self.imTitle.setObjectName("imTitle")
        #### table
        self.TauImpurityCoeffTable = QtWidgets.QTableWidget(self.tauBox)
        self.TauImpurityCoeffTable.setStyleSheet("background-color: white; border:1px solid darkGray;")
        self.TauImpurityCoeffTable.setGeometry(QtCore.QRect(187, 120, 185, 52))
        self.TauImpurityCoeffTable.setMinimumSize(QtCore.QSize(0, 0))
        self.TauImpurityCoeffTable.setFont(font_ac)
        self.TauImpurityCoeffTable.setFrameShape(QtWidgets.QFrame.Box)
        self.TauImpurityCoeffTable.setFrameShadow(QtWidgets.QFrame.Raised)
        self.TauImpurityCoeffTable.setLineWidth(1)
        self.TauImpurityCoeffTable.setMidLineWidth(1)
        self.TauImpurityCoeffTable.setAutoFillBackground(False)
        self.TauImpurityCoeffTable.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        self.TauImpurityCoeffTable.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        self.TauImpurityCoeffTable.setSelectionMode(QtWidgets.QAbstractItemView.NoSelection)
        self.TauImpurityCoeffTable.setShowGrid(True)
        self.TauImpurityCoeffTable.setGridStyle(QtCore.Qt.SolidLine)
        self.TauImpurityCoeffTable.setCornerButtonEnabled(True)
        self.TauImpurityCoeffTable.setColumnCount(3)
        self.TauImpurityCoeffTable.setRowCount(1)
        item = QtWidgets.QTableWidgetItem()
        self.TauImpurityCoeffTable.setVerticalHeaderItem(0, item)
        font_tau = QtGui.QFont()
        font_tau.setPointSize(8)
        font_tau.setBold(True)
        item = QtWidgets.QTableWidgetItem()
        item.setTextAlignment(QtCore.Qt.AlignCenter)
        item.setFont(font_tau)
        self.TauImpurityCoeffTable.setHorizontalHeaderItem(0, item)
        item = QtWidgets.QTableWidgetItem()
        item.setTextAlignment(QtCore.Qt.AlignCenter)
        item.setFont(font_tau)
        self.TauImpurityCoeffTable.setHorizontalHeaderItem(1, item)
        item = QtWidgets.QTableWidgetItem()
        item.setTextAlignment(QtCore.Qt.AlignCenter)
        item.setFont(font_tau)
        self.TauImpurityCoeffTable.setHorizontalHeaderItem(2, item)
        item = QtWidgets.QTableWidgetItem()
        item.setTextAlignment(QtCore.Qt.AlignCenter)
        self.TauImpurityCoeffTable.setItem(0, 0, item)
        item = QtWidgets.QTableWidgetItem()
        item.setTextAlignment(QtCore.Qt.AlignCenter)
        self.TauImpurityCoeffTable.setItem(0, 1, item)
        item = QtWidgets.QTableWidgetItem()
        item.setTextAlignment(QtCore.Qt.AlignCenter)
        self.TauImpurityCoeffTable.setItem(0, 2, item)
        self.TauImpurityCoeffTable.horizontalHeader().setVisible(True)
        self.TauImpurityCoeffTable.horizontalHeader().setCascadingSectionResizes(True)
        self.TauImpurityCoeffTable.horizontalHeader().setDefaultSectionSize(15)
        self.TauImpurityCoeffTable.horizontalHeader().setHighlightSections(True)
        self.TauImpurityCoeffTable.horizontalHeader().setMinimumSectionSize(15)
        self.TauImpurityCoeffTable.verticalHeader().setVisible(False)
        self.TauImpurityCoeffTable.horizontalHeader().setFixedHeight(25)
        self.TauImpurityCoeffTable.verticalHeader().setSortIndicatorShown(False)
        self.TauImpurityCoeffTable.verticalHeader().setStretchLastSection(False)
        self.TauImpurityCoeffTable.horizontalHeader().setSectionResizeMode(QtWidgets.QHeaderView.Stretch)
        self.TauImpurityCoeffTable.verticalHeader().setSectionResizeMode(QtWidgets.QHeaderView.Stretch)
        self.TauImpurityCoeffTable.setSizeAdjustPolicy(QtWidgets.QAbstractScrollArea.AdjustToContents)
        self.TauAcousticCoeffTable.itemChanged.connect(self.computetau)
        self.TauImpurityCoeffTable.setObjectName("TauImpurityCoeffTable")
        ### Matthiessen's rule
        #### title
        self.mrTitle = QtWidgets.QLabel(self.tauBox)
        self.mrTitle.setGeometry(QtCore.QRect(297, 180, 72, 20))
        font = QtGui.QFont()
        font.setPointSize(7)
        self.mrTitle.setStyleSheet("background-color: white; border:0px;")
        self.mrTitle.setFont(font)
        self.mrTitle.setObjectName("mrTitle")
        #### check box
        self.checkBoxRect = QtWidgets.QGroupBox(self.tauBox)
        self.checkBoxRect.setGeometry(QtCore.QRect(187, 197, 185, 20))
        self.checkBoxRect.setStyleSheet("background-color: white; border:2px solid {};".format(header_color))
        self.checkBoxRect.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.checkBoxRect.setObjectName("checkBoxRect")
        self.mr_checkBox_const = QtWidgets.QCheckBox(self.checkBoxRect)
        self.mr_checkBox_const.setGeometry(QtCore.QRect(4, 2, 50, 14))
        self.mr_checkBox_const.setFont(font_ac)
        self.mr_checkBox_const.setStyleSheet("border: 0px")
        self.mr_checkBox_acoust = QtWidgets.QCheckBox(self.checkBoxRect)
        self.mr_checkBox_acoust.setGeometry(QtCore.QRect(64, 2, 50, 14))
        self.mr_checkBox_acoust.setFont(font_ac)
        self.mr_checkBox_acoust.setStyleSheet("border: 0px")
        self.mr_checkBox_impur = QtWidgets.QCheckBox(self.checkBoxRect)
        self.mr_checkBox_impur.setGeometry(QtCore.QRect(124, 2, 50, 14))
        self.mr_checkBox_impur.setFont(font_ac)
        self.mr_checkBox_impur.setStyleSheet("border: 0px")
        #### gamma title
        font = QtGui.QFont()
        font.setPointSize(8)
        font.setBold(True)
        self.mr_gamma_Label = QtWidgets.QLabel(self.tauBox)
        self.mr_gamma_Label.setGeometry(QtCore.QRect(273, 219, 50, 20))
        self.mr_gamma_Label.setStyleSheet("Background-color: {}; border-top: 1px solid darkGray; border-bottom: 1px solid darkGray; border-right: 1px solid darkGray; border-left: 1px solid lightGray; color:white".format(header_color))
        self.mr_gamma_Label.setFont(font)
        self.mr_gamma_Label.setTextFormat(QtCore.Qt.PlainText)
        self.mr_gamma_Label.setAlignment(QtCore.Qt.AlignCenter)
        self.mr_gamma_Label.setObjectName("mr_gamma_Label")
        #### gamma input
        font.setBold(False)
        self.mr_gamma_Input = ClickableLineEdit(self.tauBox)
        self.mr_gamma_Input.setGeometry(QtCore.QRect(322, 219, 50, 20))
        self.mr_gamma_Input.setStyleSheet("QLineEdit {Background-color:white; border-top: 1px solid darkGray; border-bottom: 1px solid darkGray; border-right: 1px solid darkGray; border-left: 1px solid lightGray; color: gray}")
        self.mr_gamma_Input.setFont(font)
        self.mr_gamma_Input.setAlignment(QtCore.Qt.AlignCenter)
        self.mr_gamma_Input.setObjectName("muInput")
        self.mr_gamma_Input.editingFinished.connect(self.computetau)
        #### error message Matthiessen's rule
        font = QtGui.QFont()
        font.setPointSize(7)
        self.mr_error_msg = QtWidgets.QLabel(self.tauBox)
        self.mr_error_msg.setGeometry(QtCore.QRect(235, 240, 137, 12))
        self.mr_error_msg.setStyleSheet("color:red; border:0px")
        self.mr_error_msg.setFont(font)
        self.mr_error_msg.setTextFormat(QtCore.Qt.PlainText)
        self.mr_error_msg.setVisible(False)
        self.mr_error_msg.setObjectName("mr_error_msg")

        ### relaxation time plot
        self.tauplot = PlotTau(self.tauBox)
        self.tauplot.setGeometry(QtCore.QRect(1, 120, 185, 140))
        self.tauplot.setStyleSheet("background-color: white; border:1px solid {};".format(header_color))

        # import buttons
        font = QtGui.QFont()
        font.setPointSize(8)
        self.tick_icon = QtGui.QIcon()
        self.tick_icon.addPixmap(QtGui.QPixmap("./images/tick.png"), QtGui.QIcon.Active, QtGui.QIcon.Off)
        self.error_icon = QtGui.QIcon()
        self.error_icon.addPixmap(QtGui.QPixmap("./images/error.jpeg"), QtGui.QIcon.Active, QtGui.QIcon.Off)
        ## import experimental data 
        self.ImportExpButton = QtWidgets.QPushButton(self.inputBox)
        self.ImportExpButton.setGeometry(QtCore.QRect(10, 350, 120, 26))
        self.ImportExpButton.setLayoutDirection(QtCore.Qt.RightToLeft)
        self.ImportExpButton.setFont(font)
        self.ImportExpButton.setStyleSheet("QPushButton:hover{background-color: %s;}" % (header_color))
        self.ImportExpButton.clicked.connect(self.import_exp_data)
        self.ImportExpButton.setObjectName("ImportButton")
        ## clear experimental data
        self.ClearExpButton = QtWidgets.QPushButton(self.inputBox)
        self.ClearExpButton.setGeometry(QtCore.QRect(10, 375, 120, 26))
        self.ClearExpButton.setLayoutDirection(QtCore.Qt.RightToLeft)
        self.ClearExpButton.setFont(font)
        self.ClearExpButton.setStyleSheet("QPushButton:hover{background-color: %s;}" % (header_color))
        self.ClearExpButton.clicked.connect(self.clear_exp_data)
        self.ClearExpButton.setObjectName("ClearExpButton")

        # compute section
        ## progressbar
        self.progressBar = QRoundProgressBar(self.centralwidget)
        self.progressBar.setGeometry(QtCore.QRect(320, 500, 90, 90))
        self.progressBar.setBarStyle(QRoundProgressBar.StyleDonut)
        self.progressBar.setDataColors([(0., Qt_gui_color), (0.5, Qt_gui_color), (1., Qt_gui_color)])
        self.progressBar.setFixedSize(100, 100)
        self.progressBar.setNullPosition(90)
        self.progressBar.setValue(0)
        self.progressBar.setRange(0, 100)
        self.progressBar.setProperty("value", 0)
        self.progressBar.setObjectName("progressBar")
        ## compute button
        self.ComputeButton = QtWidgets.QPushButton(self.centralwidget)
        self.ComputeButton.setGeometry(330, 510, 80, 80)
        self.ComputeButton.setStyleSheet("QPushButton#ComputeButton {border-image: url(./images/grun.PNG); background-repeat: no-repeat;} QPushButton#ComputeButton:hover{border-image: url(./images/yrun.PNG); background-repeat: no-repeat;} ")
        self.ComputeButton.setObjectName("ComputeButton")
        self.ComputeButton.blockSignals(True)
        self.ComputeButton.clicked.connect(self.compute)

        self.retranslateInputUi(self.InputWindow)
        QtCore.QMetaObject.connectSlotsByName(self.InputWindow)

    def retranslateInputUi(self, InputWindow):
        _translate = QtCore.QCoreApplication.translate
        InputWindow.setWindowTitle(_translate("InputWindow", "InputWindow"))
        self.menuFile.setTitle(_translate("InputWindow", "File"))
        self.menuHelp.setTitle(_translate("InputWindow", "Help"))
        self.actionImportData.setText(_translate("InputWindow", "Import data"))
        self.actionExit.setText(_translate("InputWindow", "Exit"))
        self.actionAbout.setText(_translate("InputWindow", "About"))
        self.inputlabel.setText(_translate("InputWindow", "INPUTS"))
        self.serverLabel.setText(_translate("InputWindow", "Server status: "))
        self.ValText.setText(_translate("InputWindow", "Valence band #"))
        # item = self.ValTable.verticalHeaderItem(0)
        # item.setText(_translate("InputWindow", "New Row"))
        item = self.ValTable.horizontalHeaderItem(0)
        item.setText(_translate("InputWindow", "mx"))
        item = self.ValTable.horizontalHeaderItem(1)
        item.setText(_translate("InputWindow", "my"))
        item = self.ValTable.horizontalHeaderItem(2)
        item.setText(_translate("InputWindow", "mz"))
        item = self.ValTable.horizontalHeaderItem(3)
        item.setText(_translate("InputWindow", "E0"))
        __sortingEnabled = self.ValTable.isSortingEnabled()
        self.ValTable.setSortingEnabled(False)
        self.ValTable.setSortingEnabled(__sortingEnabled)
        self.CondText.setText(_translate("InputWindow", "Conduction band #"))
        # item = self.CondTable.verticalHeaderItem(0)
        # item.setText(_translate("InputWindow", "New Row"))
        item = self.CondTable.horizontalHeaderItem(0)
        item.setText(_translate("InputWindow", "mx"))
        item = self.CondTable.horizontalHeaderItem(1)
        item.setText(_translate("InputWindow", "my"))
        item = self.CondTable.horizontalHeaderItem(2)
        item.setText(_translate("InputWindow", "mz"))
        item = self.CondTable.horizontalHeaderItem(3)
        item.setText(_translate("InputWindow", "E0"))
        __sortingEnabled = self.CondTable.isSortingEnabled()
        self.CondTable.setSortingEnabled(False)
        self.CondTable.setSortingEnabled(__sortingEnabled)
        self.tLabel.setText(_translate("InputWindow", "T [K]"))
        self.tInput.setText(_translate("InputWindow", "[300:1000:100]"))
        self.muLabel.setText(_translate("InputWindow", "μ [eV]"))
        self.muInput.setText(_translate("InputWindow", "0.5"))
        self.tauTitle.setText(_translate("InputWindow", "𝛕 models"))
        self.comboBox.setItemText(0, _translate("MainWindow", "constant"))
        self.comboBox.setItemText(1, _translate("MainWindow", "acoustic"))
        self.comboBox.setItemText(2, _translate("MainWindow", "impurity"))
        self.comboBox.setItemText(3, _translate("MainWindow", "Matthiessen's rule"))
        self.TauplotButton.setText(_translate("MainWindow", "PLOT"))
        self.acTitle.setText(_translate("InputWindow", "Acoustic scattering"))
        self.acousBox.setItemText(0, _translate("MainWindow", " -empty- "))
        # self.acousBox.setItemText(1, _translate("MainWindow", "val_1"))
        self.imTitle.setText(_translate("InputWindow", "Impurity scattering"))
        self.mrTitle.setText(_translate("InputWindow", "Matthiessen's rule"))
        self.mr_checkBox_const.setText(_translate("InputWindow", "const"))
        self.mr_checkBox_acoust.setText(_translate("InputWindow", "acous"))
        self.mr_checkBox_impur.setText(_translate("InputWindow", "impur"))
        self.mr_gamma_Label.setText(_translate("InputWindow", "γ"))
        self.mr_gamma_Input.setText(_translate("InputWindow", "-1"))
        item = self.TauAcousticCoeffTable.horizontalHeaderItem(0)
        item.setText(_translate("InputWindow", "ε_min"))
        item = self.TauAcousticCoeffTable.item(0,0).setText('0.0')
        item = self.TauAcousticCoeffTable.horizontalHeaderItem(1)
        item.setText(_translate("InputWindow", "A_sm"))
        item = self.TauAcousticCoeffTable.item(0,1).setText('1.0')
        item = self.TauAcousticCoeffTable.horizontalHeaderItem(2)
        item.setText(_translate("InputWindow", "τm_max"))
        item = self.TauAcousticCoeffTable.item(0,2).setText('1.0')
        item = self.TauAcousticCoeffTable.horizontalHeaderItem(3)
        item.setText(_translate("InputWindow", "T₀"))
        item = self.TauAcousticCoeffTable.item(0,3).setText('50.0')
        item = self.TauAcousticCoeffTable.horizontalHeaderItem(4)
        item.setText(_translate("InputWindow", "μ_min"))
        item = self.TauAcousticCoeffTable.item(0,4).setText('2.0')
        item = self.TauAcousticCoeffTable.horizontalHeaderItem(5)
        item.setText(_translate("InputWindow", "μ_max"))
        item = self.TauAcousticCoeffTable.item(0,5).setText('2.0')
        item = self.TauImpurityCoeffTable.horizontalHeaderItem(0)
        item.setText(_translate("InputWindow", "ε_im"))
        item = self.TauImpurityCoeffTable.item(0,0).setText('0.0')
        item = self.TauImpurityCoeffTable.horizontalHeaderItem(1)
        item.setText(_translate("InputWindow", "A_im"))
        item = self.TauImpurityCoeffTable.item(0,1).setText('1.0')
        item = self.TauImpurityCoeffTable.horizontalHeaderItem(2)
        item.setText(_translate("InputWindow", "γ_im"))
        item = self.TauImpurityCoeffTable.item(0,2).setText('1.0')
        self.mr_error_msg.setText(_translate("InputWindow", "Use checkbox to select scatterings."))
        self.ComputeButton.setText(_translate("InputWindow", ""))
        self.ImportExpButton.setText(_translate("InputWindow", "IMPORT DATA"))
        self.ClearExpButton.setText(_translate("InputWindow", "CLEAR DATA"))

    # initial precalculation operations
    def check_server_status(self):
        # run first calculation to compile the code (in parallel)
        self.first_run_thread = threading.Thread(target=self.first_run, args=(), kwargs={})
        self.first_run_thread.start()
        
    # add/remove line in CondTable according to CondSpin
    @QtCore.Slot()
    def add_cond_band(self):
        self.CondTable.setRowCount(self.CondSpin.value())
        for row in range(self.CondTable.rowCount()):
            if self.CondTable.item(row, 0):
                continue
            else:
                for column in range(self.CondTable.columnCount()):
                    item = QtWidgets.QTableWidgetItem()
                    item.setTextAlignment(QtCore.Qt.AlignCenter)
                    self.CondTable.setItem(row, column, item)

    # add/remove line in ValTable according to ValSpin
    @QtCore.Slot()
    def add_val_band(self):
        self.ValTable.setRowCount(self.ValSpin.value())
        for row in range(self.ValTable.rowCount()):
            if self.ValTable.item(row, 0):
                continue
            else:
                for column in range(self.ValTable.columnCount()):
                    item = QtWidgets.QTableWidgetItem()
                    item.setTextAlignment(QtCore.Qt.AlignCenter)
                    self.ValTable.setItem(row, column, item)

    # update bands position data structure
    @QtCore.Slot()
    def update_ebandmins(self):
        self.data.ebandmins.clear()
        self.data.num_cbands = self.CondSpin.value()
        for c in range(self.data.num_cbands):
            self.data.ebandmins.append(self.CondTable.item(c, 3).text())
        self.data.num_vbands = self.ValSpin.value()
        for v in range(self.data.num_vbands):
            self.data.ebandmins.append(self.ValTable.item(v, 3).text())
        self.update_acous_box()

    @QtCore.Slot()
    def insert_new_temperature(self):
        self.tInput.setStyleSheet("QLineEdit {Background-color:white; border-top: 2px solid darkGray; border-top: 2px solid darkGray; border-bottom: 2px solid darkGray; border-left: 2px solid darkGray; border-right: 1px solid lightGray; color: black};")

    @QtCore.Slot()
    def update_temp(self):
        input_text = self.tInput.text()
        if input_text == '':
            self.tInput.setStyleSheet("QLineEdit {Background-color:white; border-top: 2px solid darkGray; border-top: 2px solid darkGray; border-bottom: 2px solid darkGray; border-left: 2px solid darkGray; border-right: 1px solid lightGray; color: gray};")
            self.tInput.setText("[300:1000:100]")

    @QtCore.Slot()
    def insert_new_mu(self):
        self.muInput.setStyleSheet("QLineEdit {Background-color:white; border-top: 2px solid darkGray; border-top: 2px solid darkGray; border-bottom: 2px solid darkGray; border-right: 2px solid darkGray; border-left: 1px solid lightGray; color: black};")

    @QtCore.Slot()
    def update_mu(self):
        input_text = self.muInput.text()
        if input_text == '':
            self.muInput.setStyleSheet("QLineEdit {Background-color:white; border-top: 2px solid darkGray; border-top: 2px solid darkGray; border-bottom: 2px solid darkGray; border-right: 2px solid darkGray; border-left: 1px solid lightGray; color: gray};")
            self.muInput.setText("0.5")


    # update acousBox with new band params
    def update_acous_box(self):
        self.acousBox.clear()
        num_cond_bands = self.CondSpin.value()
        num_val_bands = self.ValSpin.value()
        if num_cond_bands == 0 and num_val_bands == 0:
            self.acousBox.addItem(" -empty- ")
        for i in range(num_cond_bands):
            self.acousBox.addItem("cond_{}, E0 = {}".format(i+1,self.data.ebandmins[i]))
        for i in range(num_val_bands):
            self.acousBox.addItem("val_{}, E0 = {}".format(i+1,self.data.ebandmins[num_cond_bands+i]))
    

    # update all data structures
    def set_data(self):
        self.data = Data() # reset 
        self.data.num_cbands = self.CondSpin.value()
        self.data.num_vbands = self.ValSpin.value()
        self.data.num_bands = self.data.num_cbands + self.data.num_vbands
        # cond bands
        for c in range(self.data.num_cbands):
            # ctensor
            self.data.cTensors.append(self.CondTable.item(c, 0).text() + " " +
                                      self.CondTable.item(c, 1).text() + " " +
                                      self.CondTable.item(c, 2).text() + " " +
                                      "0.0" + " " + "0.0" + " " + "0.0")
            # ebandmins
            self.data.ebandmins.append(self.CondTable.item(c, 3).text())
            # mband
            self.data.mband.append("1")
            # TODO: change.This is default value
            self.data.degeneracies.append("1")
        # valence bands
        for v in range(self.data.num_vbands):
            # ctensor
            self.data.cTensors.append(self.ValTable.item(v, 0).text() + " " +
                                      self.ValTable.item(v, 1).text() + " " +
                                      self.ValTable.item(v, 2).text() + " " +
                                      "0.0" + " " + "0.0" + " " + "0.0")
            # ebandmins
            self.data.ebandmins.append(self.ValTable.item(v, 3).text())
            # mband
            self.data.mband.append("-1")
            # TODO: change.This is default value
            self.data.degeneracies.append("1")

        self.data.T_str = self.tInput.text()
        self.data.mu_str = self.muInput.text()
        
        # tau model
        self.update_tau()


    # update relaxation time data structures
    def update_tau(self):
        self.clear_tau()
        self.data.tau_model_type = self.comboBox.currentText()
        if self.data.tau_model_type == "acoustic":
            for c in range(self.TauAcousticCoeffTable.columnCount()):
                value = self.TauAcousticCoeffTable.item(0, c).text()
                if value == "":
                    self.data.tau_acoustic_coeffs.append(0.0)
                else:
                    self.data.tau_acoustic_coeffs.append(float(value))
        elif self.data.tau_model_type == "impurity":
            for c in range(self.TauImpurityCoeffTable.columnCount()):
                value = self.TauImpurityCoeffTable.item(0, c).text()
                if value == "":
                    self.data.tau_impurity_coeffs.append(0.0)
                else:
                    self.data.tau_impurity_coeffs.append(float(value))
        elif self.data.tau_model_type == "Matthiessen's rule":
            for i,tau_checkbox in enumerate([self.mr_checkBox_const,self.mr_checkBox_acoust,self.mr_checkBox_impur]):
                if tau_checkbox.isChecked():
                    models_blist = list(self.data.tau_matthiessen_models)
                    models_blist[i] = '1'
                    self.data.tau_matthiessen_models = ''.join(models_blist)
            self.data.tau_matthiessen_gamma = self.mr_gamma_Input.text()
            for c in range(self.TauAcousticCoeffTable.columnCount()):
                value = self.TauAcousticCoeffTable.item(0, c).text()
                if value == "":
                    self.data.tau_acoustic_coeffs.append(0.0)
                else:
                    self.data.tau_acoustic_coeffs.append(float(value))
            for c in range(self.TauImpurityCoeffTable.columnCount()):
                value = self.TauImpurityCoeffTable.item(0, c).text()
                if value == "":
                    self.data.tau_impurity_coeffs.append(0.0)
                else:
                    self.data.tau_impurity_coeffs.append(float(value))


    # clear relaxation time data structures
    def clear_tau(self):
        self.data.tau_model_type = "constant"
        self.data.tau_acoustic_coeffs.clear()
        self.data.tau_impurity_coeffs.clear()
        self.data.tau_matthiessen_models = "000"
        self.data.tau_matthiessen_gamma = "0.0"


    # plot relaxation time
    def publish_tau(self, data):
        print(data)
        mus = np.array(data["mu"])
        tau = np.array(data["data"])
        Ts = np.array(data["T"])
        self.tauplot.ax.cla()
        self.tauplot.plot(mus, tau, Ts)


    #### CALCULATION methods ####
    # first run to compile the code
    def first_run(self):
        while(True):
            try:
                check = requests.get('http://127.0.0.1:1200/api/check')
                check.raise_for_status() 
            except (requests.exceptions.ConnectionError, requests.exceptions.Timeout):
                continue
            except requests.exceptions.HTTPError:
                continue
            else:
                break

        url_c = 'http://127.0.0.1:1200/api/guicalc'
        headers = {'Accept': 'application/json', 'Content-Type': 'application/json'}

        # write the message
        message = {"# export all data [true/false]": False,
                    "# number of bands": 2,
                    "# Fermi level": "0.5",
                    "# temperature": '[300:500:100]',
                    "# bands masses and angles": ["1. 1. 1. 0. 0. 0.", "1. 1. 1. 0. 0. 0."],
                    "# band type": ["1","-1"],
                    "# energy extrema": ["1","-1"],
                    "# degeneracy": ["1","1"],
                    "# tau model [constant/acoustic/impurity/matthiessen]": "constant",
                    "# tau acoustic coefficients": [],
                    "# tau impurity coefficients": [],
                    "# tau matthiessen models": "000",
                    "# tau matthiessen gamma": "0.0"}

        for tensor_name in ["conductivity", "seebeck", "thermal", "concentration"]:
            message["tensor_name"] = tensor_name
            # send the request
            requests.post(url_c, json=message, headers=headers)

        # # relaxation time plot
        # self.computetau()

        # change server status signal
        self.set_greenstatus()
        # activate ComputeButton
        self.ComputeButton.blockSignals(False)
        self.TauplotButton.blockSignals(False)
        # first run completed
        self.is_first_run_completed = True


    # send a request to the server to compute the relaxation time over input ranges of temps and Fermi levels
    @QtCore.Slot()
    def computetau(self):
        # wait for the GUI to be ready
        if self.is_first_run_completed == False:
            return

        # update data structures
        self.set_data()
        self.set_redstatus()
        url_ctau = 'http://127.0.0.1:1200/api/guitaucalc'
        self.headers = {'Accept': 'application/json', 'Content-Type': 'application/json'}

        bandtype,ebandmin = 1,0

        # update acousBox
        if self.comboBox.currentText() == "acoustic" or (self.comboBox.currentText() == "Matthiessen's rule" and self.mr_checkBox_acoust):
            self.update_acous_box()
            if len(self.data.ebandmins) == 0:
                print("\033[91m[ERROR] Band structure is empty.\[\033[0m\]")
                self.set_greenstatus()
                return
            # get type of band for acoustic scattering
            acousBox = self.acousBox.currentText()
            if len(set(acousBox).intersection('cond')) == 4:
                bandtype = +1
            elif len(set(acousBox).intersection('val')) == 3:
                bandtype = -1
            else:
                raise(RuntimeError)
            # get band minimum for acoustic scattering
            ebandmin = float(acousBox.strip().split(" ")[-1])

        # write the request message
        self.message = {"# Fermi level": self.muInput.text(),
                        "# temperature": self.tInput.text(),
                        "# band type": bandtype,
                        "# energy extremum": ebandmin,
                        "# tau model [constant/acoustic/impurity/matthiessen]": self.data.tau_model_type,
                        "# tau acoustic coefficients": self.data.tau_acoustic_coeffs,
                        "# tau impurity coefficients": self.data.tau_impurity_coeffs,
                        "# tau matthiessen models": self.data.tau_matthiessen_models,
                        "# tau matthiessen gamma": self.data.tau_matthiessen_gamma}

        try:
            # send the request
            r_calc = requests.post(url_ctau, json=self.message, headers=self.headers)
            r_calc.raise_for_status()
            # check response
            if r_calc.status_code == 200:   # ok
                # check if Matthiessen had an error
                if self.mr_error_msg.isVisible():
                    self.mr_error_msg.setVisible(False)
                # plot the relaxatin time
                self.publish_tau(r_calc.json())
            # error -> clear GIU
            elif r_calc.status_code == 210 and r_calc.text == "-40":
                self.mr_error_msg.setVisible(True)
                self.clear_gui()
            elif r_calc.status_code == 210:
                self.clear_gui()
            self.set_greenstatus()

        except requests.exceptions.HTTPError as errh:
            print("Http error:", errh)
            self.set_greenstatus()
        except requests.exceptions.ConnectionError as errc:
            print ("Exception occurred. Check server.")
            self.set_greenstatus()
            raise(errc)
        except requests.exceptions.Timeout as errt:
            print("Timeout error:", errt)
            self.set_greenstatus()
        except requests.exceptions.RequestException as err:
            print("Exception occurred. Check server.", err)
            self.set_greenstatus()
            raise(err)


    # function that requests the calculation and exports the results
    @QtCore.Slot()
    def compute(self):
        self.set_redstatus()
        url_c = 'http://127.0.0.1:1200/api/guicalc'
        self.headers = {'Accept': 'application/json', 'Content-Type': 'application/json'}

        # reset progress bar
        self.progressBar.setValue(0)

        # clean the variables
        self.clear_gui()

        # update the data structure
        self.set_data()

        # write the request message
        self.message = {"# export all data [true/false]": self.data.export_all_data,
                        "# number of bands": self.data.num_bands,
                        "# Fermi level": self.data.mu_str,
                        "# temperature": self.data.T_str,
                        "# bands masses and angles": self.data.cTensors,
                        "# band type": self.data.mband,
                        "# energy extrema": self.data.ebandmins,
                        "# degeneracy": self.data.degeneracies,
                        "# tau model [constant/acoustic/impurity/matthiessen]": self.data.tau_model_type,
                        "# tau acoustic coefficients": self.data.tau_acoustic_coeffs,
                        "# tau impurity coefficients": self.data.tau_impurity_coeffs,
                        "# tau matthiessen models": self.data.tau_matthiessen_models,
                        "# tau matthiessen gamma": self.data.tau_matthiessen_gamma}

        # differently for CLI, GUI version computes all of the four tensors by default
        self.args = ["conductivity", "seebeck", "thermal", "concentration"]

        if self.is_first_run_thread_active:
            self.first_run_thread.join()
            self.is_first_run_thread_active = False

        # loop over electrical cond, Seebeck, thermal cond, carrier conc
        for idx, tensor_name in enumerate(self.args):
            self.message["tensor_name"] = tensor_name
            try:
                # send request and collect results
                r_calc = requests.post(url_c, json=self.message, headers=self.headers)
                r_calc.raise_for_status()
                # check response
                if r_calc.status_code == 200:   # ok
                    # publish results to OUTPUT window
                    self.publish_output(tensor_name, r_calc.json())
                    self.update_progress_bar(25*(idx+1))
                # error -> clear GIU
                elif r_calc.status_code == 210 and r_calc.text == "-20":
                    self.clear_gui()
                    self.ui_out.plots.figure.suptitle("", y=0.97)
                    self.ui_out.plots.ax1.spines['left'].set_visible(False)
                    self.ui_out.plots.ax1.spines['bottom'].set_visible(False)
                    self.ui_out.plots.ax1.get_xaxis().set_visible(False)
                    self.ui_out.plots.ax1.get_yaxis().set_visible(False)
                    self.ui_out.plots.ax2.spines['left'].set_visible(False)
                    self.ui_out.plots.ax2.spines['bottom'].set_visible(False)
                    self.ui_out.plots.ax2.get_xaxis().set_visible(False)
                    self.ui_out.plots.ax2.get_yaxis().set_visible(False)
                    self.ui_out.plots.figure.text(0.05, 0.85, "ERROR: ", ha="left", va="bottom", size="large", color="red", fontfamily="serif")
                    self.ui_out.plots.figure.text(0.15, 0.85, "Relaxation time functional form identically zero.", ha="left", va="bottom", size="large", fontfamily="serif")
                    self.ui_out.plots.draw()
                    self.set_greenstatus()
                    break
                elif r_calc.status_code == 210 and r_calc.text == "-30":
                    self.clear_gui()
                    self.ui_out.plots.figure.suptitle("", y=0.97)
                    self.ui_out.plots.ax1.spines['left'].set_visible(False)
                    self.ui_out.plots.ax1.spines['bottom'].set_visible(False)
                    self.ui_out.plots.ax1.get_xaxis().set_visible(False)
                    self.ui_out.plots.ax1.get_yaxis().set_visible(False)
                    self.ui_out.plots.ax2.spines['left'].set_visible(False)
                    self.ui_out.plots.ax2.spines['bottom'].set_visible(False)
                    self.ui_out.plots.ax2.get_xaxis().set_visible(False)
                    self.ui_out.plots.ax2.get_yaxis().set_visible(False)
                    self.ui_out.plots.figure.text(0.05, 0.85, "ERROR: ", ha="left", va="bottom", size="large",color="red", fontfamily="serif")
                    self.ui_out.plots.figure.text(0.15, 0.85, "Domain error in the τ function calculation.", ha="left", va="bottom", size="large", fontfamily="serif")
                    self.ui_out.plots.figure.text(0.05,0.78, "Hint:", ha="left", va="bottom", size="large",color="blue", fontfamily="serif")
                    self.ui_out.plots.figure.text(0.11,0.78, "shift the zero value of the chosen temperature or the T₀ parameter of the", ha="left", va="bottom", size="large", fontfamily="serif")
                    self.ui_out.plots.figure.text(0.11,0.71, "acoustic scattering.", ha="left", va="bottom", size="large", fontfamily="serif")
                    self.ui_out.plots.draw()
                    self.set_greenstatus()
                    break
                elif r_calc.status_code == 210 and r_calc.text == "-40":
                    self.clear_gui()

                    self.set_greenstatus()
                    break

            except requests.exceptions.HTTPError as errh:
                print("Http error:", errh)
            except requests.exceptions.ConnectionError as errc:
                print ("Exception occurred. Check server.")
                raise(errc)
            except requests.exceptions.Timeout as errt:
                print("Timeout error:", errt)
            except requests.exceptions.RequestException as err:
                print("Exception occurred. Check server.", err)
                raise(err)


    # publish output in the OUTPUT window
    def publish_output(self, tensor_name, data):

        # check if output window is visible, if not show it
        if not self.OutputWindow.isVisible():
            self.OutputWindow.show()

        # tensor has shape (6, num_mu, num_t)
        tensor = np.array(data["data"])
        norm_const = 1/3

        if tensor_name == "conductivity":

            ##### set slider for T and mu #####
            self.T = np.array(data["T"])
            self.mus = np.array(data["mu"])
            if self.T.shape == ():
                self.T = np.array([self.T])
            if self.mus.shape == ():
                self.mus = np.array([self.mus])
            self.num_t = self.T.size
            self.num_mu = self.mus.size
            if not self.out_all_data.isallocated():
                self.out_all_data.setTmu(self.num_mu, self.num_t)
            if self.T.size < 2:
                # set slider according to user inputs
                self.ui_out.TSlider.setMinimum(self.T[0])
                self.ui_out.TSlider.setMaximum(self.T[0])
                self.ui_out.TSlider.setSingleStep(0)
            else:
                # set slider according to user inputs
                self.ui_out.TSlider.setMinimum(0)
                self.ui_out.TSlider.setMaximum(len(self.T)-1)
                self.T_step = self.T[1] - self.T[0]
                self.ui_out.TSlider.setSingleStep(1)
                self.ui_out.TSlider.setTickInterval(1)
                self.ui_out.TSlider.setTickPosition(QtWidgets.QSlider.TicksBelow)
            if self.mus.size < 2:
                # set slider according to user inputs
                self.ui_out.muSlider.setMinimum(self.mus[0])
                self.ui_out.muSlider.setMaximum(self.mus[0])
                self.ui_out.muSlider.setSingleStep(0)
                self.muStepConv = 0.0
            else:
                # set slider according to user inputs
                self.mu_step = self.data.mu_str.split(":")[-1]
                self.ui_out.muSlider.setMinimum(0)
                self.ui_out.muSlider.setMaximum(len(self.mus)-1)
                self.ui_out.muSlider.setSingleStep(1)
                self.ui_out.muSlider.setTickInterval(1)
                self.ui_out.muSlider.setTickPosition(QtWidgets.QSlider.TicksBelow)

            # set T mu values in TVal and muVal
            self.ui_out.TVal.setText(str(int(self.T[0])))
            self.ui_out.muVal.setText(str(self.mus[0]))
            # activate sliders
            self.ui_out.TVal.textChanged.connect(self.ui_out.TValueChanged)
            self.ui_out.muVal.textChanged.connect(self.ui_out.muValueChanged)
            self.ui_out.TSlider.valueChanged.connect(self.ui_out.TSliderChanged)
            self.ui_out.muSlider.valueChanged.connect(self.ui_out.muSliderChanged)
            
        ##### compute the trace for each tensor
            self.out_all_data.setCond(tensor)
            trace_tensor = np.empty((self.num_mu, self.num_t))
            for t in range(self.num_t):
                for m in range(self.num_mu):
                    trace_tensor[m, t] = np.multiply((tensor[0, m, t]+tensor[1, m, t]+tensor[2, m, t]), norm_const)
        elif tensor_name == "seebeck":
            self.out_all_data.setSeebeck(tensor)
            trace_tensor = np.empty((self.mus.size, self.T.size))
            for t in range(self.T.size):
                for m in range(self.mus.size):
                    trace_tensor[m, t] = np.multiply((tensor[0, m, t] + tensor[1, m, t] + tensor[2, m, t]), norm_const)
        elif tensor_name == "thermal":
            self.out_all_data.setThermal(tensor)
            trace_tensor = np.empty((self.mus.size, self.T.size))
            for t in range(self.T.size):
                for m in range(self.mus.size):
                    trace_tensor[m, t] = np.multiply((tensor[0, m, t]+tensor[1, m, t]+tensor[2, m, t]), norm_const)
        elif tensor_name == "concentration":
            self.out_all_data.setConc(tensor)
            trace_tensor = tensor[0, :, :]

        # plot the results
        self.ui_out.plots.plot(tensor_name, self.T, trace_tensor, self.mus, self.data.tau_model_type, self.exp_data)
        
        self.out_trace_data.data.append(trace_tensor)
        self.out_trace_data.label.append(tensor_name)

        # after the last tensor is plotted -> green light
        if tensor_name == "concentration":
            self.set_greenstatus()


    # publish each tensor in the outputTable according to the sliders (or TVal and muVal)
    def publish_tensor(self):
        mu = float(self.ui_out.muVal.text())
        t = float(self.ui_out.TVal.text())
        if self.T.shape == ():
            self.T = np.array([self.T])
        if self.mus.shape == ():
            self.mus = np.array([self.mus])
        if t not in self.T or mu not in self.mus:
            return
        if (self.ui_out.TVal.text() == "") or (self.ui_out.muVal.text() == "") \
                or (mu < self.mus[0]) or (mu > self.mus[-1]) \
                or (t < self.T[0]) or (t > self.T[-1]):
            return
        else:
            mu_idx = np.where(self.mus == mu)[0][0]
            t_idx = np.where(self.T == t)[0][0]

            tensor_idx = 0
            if self.out_all_data.conductivity is not None and self.out_all_data.conductivity.size != 0:
                self.update_row_outputTable(0, self.out_all_data.conductivity[:, mu_idx, t_idx])
                y = self.out_trace_data.data[tensor_idx][mu_idx][t_idx]
                if self.ui_out.plots.point1 is not None:
                    self.ui_out.plots.point1.remove()
                self.ui_out.plots.point1, = self.ui_out.plots.ax1.plot(t, y, marker='.', color="#1f77b4" if self.mus.size == 1 else 'dimgray', zorder=10)
                tensor_idx += 1
            if self.out_all_data.seebeck is not None and self.out_all_data.seebeck.size != 0:
                self.update_row_outputTable(1, np.multiply(self.out_all_data.seebeck[:, mu_idx, t_idx],1e6))
                y = self.out_trace_data.data[tensor_idx][mu_idx][t_idx]
                if self.ui_out.plots.point2 is not None:
                    self.ui_out.plots.point2.remove()
                self.ui_out.plots.point2, = self.ui_out.plots.ax2.plot(t, np.multiply(y,1e6), marker='.', color="orange" if self.mus.size == 1 else 'dimgray', zorder=10)
                tensor_idx += 1
            if self.out_all_data.thermal is not None and self.out_all_data.thermal.size != 0:
                self.update_row_outputTable(2, self.out_all_data.thermal[:, mu_idx, t_idx])
                y = self.out_trace_data.data[tensor_idx][mu_idx][t_idx]
                if self.ui_out.plots.point3 is not None:
                    self.ui_out.plots.point3.remove()
                self.ui_out.plots.point3, = self.ui_out.plots.ax3.plot(t, y, marker='.', color="red" if self.mus.size == 1 else 'dimgray', markersize=3, zorder=10)
                tensor_idx += 1
            if self.out_all_data.concentration is not None and self.out_all_data.concentration.size != 0:
                self.update_n_outputTable(self.out_all_data.concentration[:, mu_idx, t_idx])
                y = self.out_trace_data.data[tensor_idx][mu_idx][t_idx]
                if self.ui_out.plots.point4 is not None:
                    self.ui_out.plots.point4.remove()
                self.ui_out.plots.point4, = self.ui_out.plots.ax4.plot(t, y, marker='.', color="limegreen" if self.mus.size == 1 else 'dimgray', zorder=10)
        self.ui_out.plots.draw()


    def update_progress_bar(self, perc):
        self.progressBar.setValue(perc)


    # publish conductivities and Seebeck on outputTable
    def update_row_outputTable(self, row, data):
        self.ui_out.outputTable.item(row, 0).setText(format(data[0], '.5e'))
        self.ui_out.outputTable.item(row, 1).setText(format(data[1], '.5e'))
        self.ui_out.outputTable.item(row, 2).setText(format(data[2], '.5e'))
        self.ui_out.outputTable.item(row, 3).setText(format(data[3], '.5e'))
        self.ui_out.outputTable.item(row, 4).setText(format(data[4], '.5e'))
        self.ui_out.outputTable.item(row, 5).setText(format(data[5], '.5e'))

    # carrier concentration is a scalar 
    def update_n_outputTable(self, data):
        self.ui_out.outputTable.item(3, 0).setText(format(data[0], '.5e'))
        self.ui_out.outputTable.item(3, 1).setText('-')
        self.ui_out.outputTable.item(3, 2).setText('-')
        self.ui_out.outputTable.item(3, 3).setText('-')
        self.ui_out.outputTable.item(3, 4).setText('-')
        self.ui_out.outputTable.item(3, 5).setText('-')


    @QtCore.Slot()
    def import_exp_data(self):
        fileName, selectedFilter = QtWidgets.QFileDialog.getOpenFileName(self.InputWindow, 'Import File', str(os.getcwd()),
                                                             "CSV (Comma delimited) (*.csv);; Excel Workbook (*.xlsx)")
        try:
            if selectedFilter == "CSV (Comma delimited) (*.csv)":
                self.exp_data.set_exp_data(pd.read_csv(fileName, header=None))
            elif selectedFilter == "Excel Workbook (*.xlsx)":
                self.exp_data.set_exp_data(pd.read_excel(fileName, header=None))
            self.ImportExpButton.setIcon(self.tick_icon)
        except ValueError:
            self.ImportExpButton.setIcon(self.error_icon)
        except FileNotFoundError:
            pass
        self.ClearExpButton.setIcon(QtGui.QIcon())


    @QtCore.Slot()
    def clear_exp_data(self):
        self.exp_data.clear()
        self.ClearExpButton.setIcon(self.tick_icon)
        self.ImportExpButton.setIcon(QtGui.QIcon())


    def export_data(self, filename):
        hor_header = ['μ']
        for i in range(self.num_t):
            hor_header.append('T_'+str(i+1))
        # list not empty
        if self.out_trace_data.data:
            tensor_num = len(self.out_trace_data.data)
            mu_size, t_size = self.out_trace_data.data[0].shape
            data_to_write = np.empty((mu_size*tensor_num, len(hor_header)))
            for t in range(tensor_num):
                tensor = np.insert(self.out_trace_data.data[t], 0, self.mus, axis=1)
                data_to_write[mu_size*t:mu_size*(t+1), :] = tensor

            lbs = list(itertools.chain.from_iterable(itertools.repeat(x, mu_size) for x in self.out_trace_data.label))
            pd.DataFrame(np.asarray(data_to_write), columns=hor_header, index=lbs).to_csv(filename)


    # clear all datastructures and plots
    def clear_gui(self):
        self.data = Data()
        self.out_trace_data.clear()
        self.out_all_data.clear()
        self.ui_out.plots.ax1.cla(); self.ui_out.plots.ax2.cla(); self.ui_out.plots.ax3.cla(); self.ui_out.plots.ax4.cla();
        # space for error message. TODO
        for text in self.ui_out.plots.figure.texts:
            text.set_visible(False)
        self.ui_out.plots.ax1.spines['left'].set_visible(True)
        self.ui_out.plots.ax1.spines['bottom'].set_visible(True)
        self.ui_out.plots.ax1.get_yaxis().set_visible(True)
        self.ui_out.plots.ax1.get_yaxis().set_visible(True)
        self.ui_out.plots.ax2.spines['left'].set_visible(True)
        self.ui_out.plots.ax2.spines['bottom'].set_visible(True)
        self.ui_out.plots.ax2.get_xaxis().set_visible(True)
        self.ui_out.plots.ax2.get_yaxis().set_visible(True)
        self.ui_out.plots.draw()
        for r in range(self.ui_out.outputTable.rowCount()):
            for c in range(self.ui_out.outputTable.columnCount()):
                item = QtWidgets.QTableWidgetItem()
                item.setTextAlignment(QtCore.Qt.AlignCenter)
                self.ui_out.outputTable.setItem(r, c, item)


    def exit(self):
        self.ui_out.OutputWindow.close()
        self.InputWindow.close()


    def set_redstatus(self):
        # down signal to on
        self.redSignal.setVisible(True)
        # up signal to off
        self.greenSignal.setVisible(False)
        QtWidgets.QApplication.processEvents()


    def set_greenstatus(self):
        # down signal to off
        self.redSignal.setVisible(False)
        # activate up signal
        self.greenSignal.setVisible(True)
        QtWidgets.QApplication.processEvents()



# output window
class UiOutputWindow(object):
    def __init__(self, parent):
        self.parent = parent

    def setupUi(self, OutputWindow):
        self.OutputWindow = OutputWindow
        self.OutputWindow.setObjectName("OutputWindow")
        self.OutputWindow.resize(730, 750)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.OutputWindow.sizePolicy().hasHeightForWidth())
        self.OutputWindow.setSizePolicy(sizePolicy)
        self.OutputWindow.setAutoFillBackground(False)
        self.OutputWindow.setWindowIcon(self.parent.appIcon)
        self.plots = PlotsCanvas(width=20, height=50)
        self.centralwidget = QtWidgets.QWidget(self.OutputWindow)
        self.centralwidget.setObjectName("centralwidget")
        self.OutputWindow.setCentralWidget(self.centralwidget)

        # menubar
        self.menubar = QtWidgets.QMenuBar(self.OutputWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 819, 20))
        self.menubar.setObjectName("menubar")
        self.menuFile = QtWidgets.QMenu(self.menubar)
        self.menuFile.setObjectName("menuFile")
        self.menuHelp = QtWidgets.QMenu(self.menubar)
        self.menuHelp.setObjectName("menuHelp")
        self.OutputWindow.setMenuBar(self.menubar)
        self.statusbar = QtWidgets.QStatusBar(self.OutputWindow)
        self.statusbar.setObjectName("statusbar")
        self.OutputWindow.setStatusBar(self.statusbar)
        self.actionSave_plots = QtWidgets.QAction(self.OutputWindow)
        self.actionSave_plots.setObjectName("actionSave_plots")
        self.actionSave_plots.triggered.connect(self.create_saveplot_dialog)
        self.actionSave_plots.setShortcut("Ctrl+S")
        self.actionExport_data = QtWidgets.QAction(self.OutputWindow)
        self.actionExport_data.setObjectName("actionExport_data")
        self.actionExport_data.setStatusTip('Save File')
        self.actionExport_data.triggered.connect(self.save_data)
        self.actionExit = QtWidgets.QAction(self.OutputWindow)
        self.actionExit.setObjectName("actionExit")
        self.actionExit.triggered.connect(self.OutputWindow.close)
        self.actionAbout = QtWidgets.QAction(self.OutputWindow)
        self.actionAbout.setObjectName("actionAbout")
        self.menubar.addAction(self.menuFile.menuAction())
        self.menubar.addAction(self.menuHelp.menuAction())
        self.menuFile.addAction(self.actionSave_plots)
        self.menuFile.addAction(self.actionExport_data)
        self.menuFile.addAction(self.actionExit)
        self.menuHelp.addAction(self.actionAbout)

        # layout grid to handle outputs
        self.frameLeft = QtWidgets.QFrame(self.centralwidget)
        self.frameLeft.setMinimumSize(QtCore.QSize(0, 70))
        self.frameLeft.setMaximumSize(QtCore.QSize(16777215, 70))
        self.frameLeft.setFrameShape(QtWidgets.QFrame.NoFrame)
        self.frameLeft.setFrameShadow(QtWidgets.QFrame.Plain)
        self.frameLeft.setObjectName("frameLeft")
        self.frameTmu = QtWidgets.QFrame(self.centralwidget)
        self.frameTmu.setMinimumSize(QtCore.QSize(530, 70))
        self.frameTmu.setMaximumSize(QtCore.QSize(16777215, 70))
        self.frameTmu.setFrameShape(QtWidgets.QFrame.NoFrame)
        self.frameTmu.setFrameShadow(QtWidgets.QFrame.Plain)
        self.frameTmu.setObjectName("frameTmu")
        self.gridLayoutPlots = QtWidgets.QGridLayout()
        self.gridLayoutPlots.setObjectName("gridLayoutPlots")
        self.TmuLayout = QtWidgets.QHBoxLayout()
        self.TmuLayout.addWidget(self.frameLeft)
        self.TmuLayout.setObjectName("TmuLayout")
        self.TmuLayout.setContentsMargins(0, 20, 0, 0)
        self.frameRight = QtWidgets.QFrame(self.centralwidget)
        self.frameRight.setMaximumSize(QtCore.QSize(16777215, 70))
        self.frameRight.setFrameShape(QtWidgets.QFrame.NoFrame)
        self.frameRight.setFrameShadow(QtWidgets.QFrame.Plain)
        self.frameRight.setObjectName("frameRight")
        self.TmuLayout.addWidget(self.frameRight)
        self.gridLayoutTensors = QtWidgets.QGridLayout()
        self.gridLayoutTensors.setObjectName("gridLayoutTensors")
        self.gridLayoutTensors.setContentsMargins(50, 0, 50, 50)
        self.windowLayout = QtWidgets.QGridLayout(self.centralwidget)
        self.windowLayout.setObjectName("gridLayout_2")
        self.windowLayout.addLayout(self.gridLayoutPlots, 0, 0, 1, 1)
        self.windowLayout.addLayout(self.TmuLayout, 1, 0, 1, 1)
        self.windowLayout.addLayout(self.gridLayoutTensors, 6, 0, 1, 1)

        # title
        font_title = QtGui.QFont()
        font_title.setPointSize(15)
        font_title.setBold(True)
        font_title.setWeight(75)
        self.outputlabel = QtWidgets.QLabel(self.centralwidget)
        self.outputlabel.setFont(font_title)
        self.outputlabel.setTextFormat(QtCore.Qt.PlainText)
        self.outputlabel.setAlignment(QtCore.Qt.AlignCenter)
        self.outputlabel.setMaximumSize(QtCore.QSize(16777215, 30))
        self.outputlabel.setIndent(0)
        self.outputlabel.setObjectName("outputlabel")
        self.gridLayoutPlots.addWidget(self.outputlabel, 0, 0, 1, 1)

        # Plot table
        self.gridLayoutPlots.addWidget(self.plots, 1, 0, 1, 1)

        ## Thermodynamics params
        font = QtGui.QFont()
        font.setPointSize(9)
        TmuLabelSizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Maximum)
        TmuLabelSizePolicy.setHorizontalStretch(0)
        TmuLabelSizePolicy.setVerticalStretch(0)
        ### temperature input namebox
        self.TSym = QtWidgets.QLabel(self.frameTmu)
        self.TSym.setObjectName("TSym")
        self.TSym.setGeometry(QtCore.QRect(30, 10, 30, 16))
        self.TSym.setMaximumSize(QtCore.QSize(16777215, 20))
        self.TSym.setFont(font)
        sizePolicy.setHeightForWidth(self.TSym.sizePolicy().hasHeightForWidth())
        self.TSym.setSizePolicy(sizePolicy)
        ### temperature input 
        self.TVal = QtWidgets.QLineEdit(self.frameTmu)
        self.TVal.setObjectName("TmuVal")
        self.TVal.setGeometry(QtCore.QRect(70, 10, 71, 16))
        self.TVal.setFont(font)
        self.TVal.setAlignment(QtCore.Qt.AlignRight)
        self.TVal.returnPressed.connect(self.TValueChangedReturn)
        ### temperature slider
        self.TSlider = QtWidgets.QSlider(self.frameTmu)
        self.TSlider.setObjectName("TSlider")
        self.TSlider.setGeometry(QtCore.QRect(10, 40, 170, 20))
        self.TSlider.setMaximumSize(QtCore.QSize(16777215, 20))
        self.TSlider.setOrientation(QtCore.Qt.Horizontal)
        self.TSlider.setStyleSheet("QSlider::handle:horizontal {background-color: %s;}" % (header_color))
        ### Fermi level input namebox
        self.muSym = QtWidgets.QLabel(self.frameTmu)
        self.muSym.setObjectName("muSym")
        self.muSym.setGeometry(QtCore.QRect(220, 10, 40, 16))
        self.muSym.setMaximumSize(QtCore.QSize(16777215, 20))
        self.muSym.setFont(font)
        sizePolicy.setHeightForWidth(self.muSym.sizePolicy().hasHeightForWidth())
        self.muSym.setSizePolicy(sizePolicy)
        ### Fermi level input 
        self.muVal = QtWidgets.QLineEdit(self.frameTmu)
        self.muVal.setObjectName("muVal")
        self.muVal.setGeometry(QtCore.QRect(270, 10, 71, 16))
        self.muVal.setFont(font)
        self.muVal.setAlignment(QtCore.Qt.AlignRight)
        self.muVal.returnPressed.connect(self.muValueChangedReturn)
        ### Fermi level slider
        self.muSlider = QtWidgets.QSlider(self.frameTmu)
        self.muSlider.setObjectName("muSlider")
        self.muSlider.setGeometry(QtCore.QRect(200, 40, 170, 20))
        self.muSlider.setMaximumSize(QtCore.QSize(16777215, 20))
        self.muSlider.setOrientation(QtCore.Qt.Horizontal)
        self.muSlider.setStyleSheet("QSlider::handle:horizontal {background-color: %s;}" % (header_color))
        self.TmuLayout.addWidget(self.frameTmu)

        ## outputs table
        self.outputTable = QtWidgets.QTableWidget(self.centralwidget)
        self.outputTable.setObjectName("outputTable")
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.outputTable.sizePolicy().hasHeightForWidth())
        sizePolicy.setWidthForHeight(self.outputTable.sizePolicy().hasWidthForHeight())
        self.outputTable.setSizePolicy(sizePolicy)
        self.outputTable.setMinimumSize(QtCore.QSize(620, 117))
        self.outputTable.setMaximumSize(QtCore.QSize(620, 117))
        self.outputTable.horizontalHeader().setDefaultSectionSize(100)
        self.outputTable.verticalHeader().setDefaultSectionSize(15)
        self.outputTable.verticalHeader().setMinimumSectionSize(15)
        self.outputTable.horizontalHeader().setFixedHeight(22)
        self.outputTable.setSizeAdjustPolicy(QtWidgets.QAbstractScrollArea.AdjustToContents)
        self.outputTable.setShowGrid(True)
        self.outputTable.setObjectName("outputTable")
        otfont = QtGui.QFont()
        otfont.setBold(False)
        otfont.setPointSize(9)
        self.outputTable.setFont(otfont)
        self.outputTable.setColumnCount(6)
        self.outputTable.setRowCount(4)
        item = QtWidgets.QTableWidgetItem()
        item.setTextAlignment(QtCore.Qt.AlignCenter)
        self.outputTable.setVerticalHeaderItem(0, item)
        item = QtWidgets.QTableWidgetItem()
        item.setTextAlignment(QtCore.Qt.AlignCenter)
        self.outputTable.setVerticalHeaderItem(1, item)
        item = QtWidgets.QTableWidgetItem()
        item.setTextAlignment(QtCore.Qt.AlignCenter)
        self.outputTable.setVerticalHeaderItem(2, item)
        item = QtWidgets.QTableWidgetItem()
        item.setTextAlignment(QtCore.Qt.AlignCenter)
        self.outputTable.setVerticalHeaderItem(3, item)
        item = QtWidgets.QTableWidgetItem()
        item.setTextAlignment(QtCore.Qt.AlignCenter)
        self.outputTable.setHorizontalHeaderItem(0, item)
        item = QtWidgets.QTableWidgetItem()
        item.setTextAlignment(QtCore.Qt.AlignCenter)
        self.outputTable.setHorizontalHeaderItem(1, item)
        item = QtWidgets.QTableWidgetItem()
        item.setTextAlignment(QtCore.Qt.AlignCenter)
        self.outputTable.setHorizontalHeaderItem(2, item)
        item = QtWidgets.QTableWidgetItem()
        item.setTextAlignment(QtCore.Qt.AlignCenter)
        self.outputTable.setHorizontalHeaderItem(3, item)
        item = QtWidgets.QTableWidgetItem()
        item.setTextAlignment(QtCore.Qt.AlignCenter)
        self.outputTable.setHorizontalHeaderItem(4, item)
        item = QtWidgets.QTableWidgetItem()
        item.setTextAlignment(QtCore.Qt.AlignCenter)
        self.outputTable.setHorizontalHeaderItem(5, item)
        self.outputTable.horizontalHeader().setSectionResizeMode(QtWidgets.QHeaderView.Stretch)
        self.outputTable.verticalHeader().setSectionResizeMode(QtWidgets.QHeaderView.Stretch)
        self.gridLayoutTensors.addWidget(self.outputTable, 0, 0, 1, 1)

        self.retranslateOutputUi(self.OutputWindow)
        QtCore.QMetaObject.connectSlotsByName(self.OutputWindow)

    def retranslateOutputUi(self, OutputWindow):
        _translate = QtCore.QCoreApplication.translate
        OutputWindow.setWindowTitle(_translate("OutputWindow", "OutputWindow"))
        self.outputlabel.setText(_translate("OutputWindow", "OUTPUTS"))
        self.TSym.setText(_translate("OutputWindow", "T [K]:"))
        self.TVal.setText(_translate("OutputWindow", "0"))
        self.muSym.setText(_translate("OutputWindow", "μ [eV]:")) # TODO
        self.muVal.setText(_translate("OutputWindow", "0.0"))
        item = self.outputTable.verticalHeaderItem(0)
        item.setText(_translate("OutputWindow", "σ"))
        item = self.outputTable.verticalHeaderItem(1)
        item.setText(_translate("OutputWindow", "S"))
        item = self.outputTable.verticalHeaderItem(2)
        item.setText(_translate("OutputWindow", "κ" + u"\u2091"))
        item = self.outputTable.verticalHeaderItem(3)
        item.setText(_translate("OutputWindow", "n"))
        item = self.outputTable.horizontalHeaderItem(0)
        item.setText(_translate("OutputWindow", "11"))
        item = self.outputTable.horizontalHeaderItem(1)
        item.setText(_translate("OutputWindow", "22"))
        item = self.outputTable.horizontalHeaderItem(2)
        item.setText(_translate("OutputWindow", "33"))
        item = self.outputTable.horizontalHeaderItem(3)
        item.setText(_translate("OutputWindow", "12"))
        item = self.outputTable.horizontalHeaderItem(4)
        item.setText(_translate("OutputWindow", "13"))
        item = self.outputTable.horizontalHeaderItem(5)
        item.setText(_translate("OutputWindow", "23"))
        self.menuFile.setTitle(_translate("OutputWindow", "File"))
        self.menuHelp.setTitle(_translate("OutputWindow", "Help"))
        self.actionSave_plots.setText(_translate("OutputWindow", "Save plots"))
        self.actionExport_data.setText(_translate("OutputWindow", "Export data"))
        self.actionExit.setText(_translate("InputWindow", "Exit"))
        self.actionAbout.setText(_translate("OutputWindow", "About"))


    @QtCore.Slot()
    def create_saveplot_dialog(self):
        self.SavePlotDialog = QtWidgets.QDialog()
        self.ui_savePlot = UiSavePlotDialog(self)
        self.ui_savePlot.setupUi(self.SavePlotDialog)
        self.SavePlotDialog.show()


    # save calculations
    @QtCore.Slot()
    def save_data(self):
        filename, _ = QtWidgets.QFileDialog.getSaveFileName(self.OutputWindow, 'Save File', str(os.getcwd()),
                                                             "CSV Files (*.csv)")
        if filename[-3:] != "csv":
            filename = filename + ".csv"
        self.parent.export_data(filename)


    # update outable when T slider changes
    @QtCore.Slot()
    def TValueChanged(self):
        try:
            if float(self.TVal.text()) not in self.parent.T:
                return
            else:
                self.parent.publish_tensor()
        except ValueError:
            print("Wrong input number.")


    # update outable when Fermi level slider changes
    @QtCore.Slot()
    def muValueChanged(self):
        try:
            if float(self.muVal.text()) not in self.parent.mus:
                return
            else:
                self.parent.publish_tensor()
        except ValueError:
            print("Wrong input number.")


    # update outable when T label changes
    @QtCore.Slot()
    def TValueChangedReturn(self):
        if float(self.TVal.text()) not in self.parent.T:
            return
        if self.parent.T.size > 1:
            value = int(self.TVal.text()-self.parent.T[0])/self.parent.T_step
            self.TSlider.setValue(int(value))
        self.parent.publish_tensor()


    # update outable when Fermi level label changes
    @QtCore.Slot()
    def muValueChangedReturn(self):
        if float(self.muVal.text()) not in self.parent.mus:
            return
        if self.parent.mus.size > 1:
            rounding_digits = decimal_digits(float(self.parent.mu_step))
            value = round((float(self.muVal.text()) - self.parent.mus[0]) / float(self.parent.mu_step), rounding_digits)
            self.muSlider.setValue(int(value))
        self.parent.publish_tensor()

    # update T label when slider changes
    @QtCore.Slot()
    def TSliderChanged(self):
        value = int(self.TSlider.value()*self.parent.T_step + self.parent.T[0])
        self.TVal.setText(str(value))

    # udapte Fermi level labels when slider changes
    @QtCore.Slot()
    def muSliderChanged(self):
        rounding_digits = decimal_digits(self.parent.mu_step)
        value = round(self.muSlider.value()*float(self.parent.mu_step) + self.parent.mus[0], rounding_digits)
        self.muVal.setText(str(value))


    @QtCore.Slot()
    def set_last_tick(self):
        self.TSlider.last_tick_pos = self.TSlider.value()


    @QtCore.Slot()
    def check_values(self):
        self.TSlider.check_values(self.TSlider.value())


# save plot dialog
class UiSavePlotDialog(object):
    def __init__(self, parent):
        self.parent = parent

    def setupUi(self, saveDialog):
        self.saveDialog = saveDialog
        self.saveDialog.setObjectName("Dialog")
        self.saveDialog.resize(536, 164)
        self.gridLayoutSave = QtWidgets.QGridLayout(self.saveDialog)
        self.gridLayoutSave.setObjectName("gridLayoutSave")

        # set font
        font = QtGui.QFont()
        font.setPointSize(9)
        # set common policy
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        # export path label
        self.pathLabel = QtWidgets.QLabel(self.saveDialog)
        self.pathLabel.setObjectName("pathLabel")
        self.pathLabel.setFont(font)
        self.pathLabel.setMaximumSize(QtCore.QSize(75, 16777215))
        sizePolicy.setHeightForWidth(self.pathLabel.sizePolicy().hasHeightForWidth())
        self.pathLabel.setSizePolicy(sizePolicy)
        self.pathLabel.setFrameShape(QtWidgets.QFrame.Box)
        self.pathLabel.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.gridLayoutSave.addWidget(self.pathLabel, 0, 0, 1, 1)
        # export path input box
        self.pathInput = QtWidgets.QLineEdit(self.saveDialog)
        self.pathInput.setFont(font)
        self.pathInput.setObjectName("lineEdit_3")
        self.gridLayoutSave.addWidget(self.pathInput, 0, 1, 1, 1)
        # browse button
        self.browseButton = QtWidgets.QPushButton(self.saveDialog)
        self.browseButton.setObjectName("BrowseButton")
        self.browseButton.setFont(font)
        self.browseButton.setMaximumSize(QtCore.QSize(50, 16777215))
        self.browseButton.setLayoutDirection(QtCore.Qt.LeftToRight)
        self.browseButton.clicked.connect(self.browse)
        self.gridLayoutSave.addWidget(self.browseButton, 0, 2, 1, 1)
        # filename label
        self.filenameLabel = QtWidgets.QLabel(self.saveDialog)
        self.filenameLabel.setObjectName("filenameLabel")
        self.filenameLabel.setFont(font)
        self.filenameLabel.setMaximumSize(QtCore.QSize(75, 16777215))
        sizePolicy.setHeightForWidth(self.filenameLabel.sizePolicy().hasHeightForWidth())
        self.filenameLabel.setSizePolicy(sizePolicy)
        self.filenameLabel.setFrameShape(QtWidgets.QFrame.Box)
        self.filenameLabel.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.gridLayoutSave.addWidget(self.filenameLabel, 2, 0, 1, 1)
        # filename inputbox
        self.filenameInput = QtWidgets.QLineEdit(self.saveDialog)
        self.filenameInput.setFont(font)
        self.filenameInput.setObjectName("lineEdit_2")
        self.gridLayoutSave.addWidget(self.filenameInput, 2, 1, 1, 1)
        # format
        self.formatBox = QtWidgets.QComboBox(self.saveDialog)
        self.formatBox.setObjectName("formatBox")
        self.formatBox.setFont(font)
        self.formatBox.addItem("jpg")
        self.formatBox.addItem("png")
        self.formatBox.addItem("pdf")
        self.formatBox.setMaximumSize(QtCore.QSize(50, 16777215))
        self.gridLayoutSave.addWidget(self.formatBox, 2, 2, 1, 1)
        # dpi label
        self.dpiLabel = QtWidgets.QLabel(self.saveDialog)
        self.dpiLabel.setObjectName("dpiLabel")
        self.dpiLabel.setFont(font)
        self.dpiLabel.setMaximumSize(QtCore.QSize(75, 16777215))
        sizePolicy.setHeightForWidth(self.dpiLabel.sizePolicy().hasHeightForWidth())
        self.dpiLabel.setSizePolicy(sizePolicy)
        self.dpiLabel.setFrameShape(QtWidgets.QFrame.Box)
        self.dpiLabel.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.gridLayoutSave.addWidget(self.dpiLabel, 3, 0, 1, 1)
        # dpi input box
        self.dpiInput = ClickableLineEdit(self.saveDialog)
        self.dpiInput.setObjectName("lineEdit")
        self.dpiInput.setStyleSheet("color: gray")
        self.dpiInput.setFont(font)
        self.dpiInput.clicked.connect(self.insert_new_dpi)
        self.dpiInput.editingFinished.connect(self.update_dpi)
        self.gridLayoutSave.addWidget(self.dpiInput, 3, 1, 1, 1)
        # save button
        self.saveButton = QtWidgets.QPushButton(self.saveDialog)
        font = QtGui.QFont()
        font.setPointSize(9)
        font.setBold(True)
        self.saveButton.setFont(font)
        self.saveButton.setLayoutDirection(QtCore.Qt.LeftToRight)
        self.saveButton.setObjectName("SaveButton")
        self.saveButton.clicked.connect(self.save_plots)
        self.gridLayoutSave.addWidget(self.saveButton, 4, 0, 1, 3)

        self.retranslateUi(self.saveDialog)
        QtCore.QMetaObject.connectSlotsByName(self.saveDialog)

    def retranslateUi(self, Dialog):
        _translate = QtCore.QCoreApplication.translate
        Dialog.setWindowTitle(_translate("Dialog", "Save plots"))
        self.pathLabel.setText(_translate("Dialog", "Path"))
        self.filenameLabel.setText(_translate("Dialog", "Filename"))
        self.browseButton.setText(_translate("Dialog", "..."))
        self.dpiLabel.setText(_translate("Dialog", "dpi"))
        self.dpiInput.setText(_translate("Dialog", "300"))
        self.saveButton.setText(_translate("Dialog", "Save"))

    @QtCore.Slot()
    def save_plots(self):
        self.path = self.pathInput.text()
        self.filename = self.filenameInput.text()
        self.dpi = float(self.dpiInput.text())
        self.format = self.formatBox.currentText()
        if (self.filename[-3:] != "jpg") or (self.filename[-3:] != "png") or (self.filename[-3:] != "pdf"):
            self.filename = self.filename + "." + self.format
        self.fullpath = self.path + "/" + self.filename
        self.parent.plots.save(self.fullpath, self.dpi)
        self.saveDialog.close()

    @QtCore.Slot()
    def browse(self):
        path = QtWidgets.QFileDialog.getExistingDirectory(self.parent.OutputWindow, 'Save File', str(os.getcwd()))
        self.pathInput.setText(path)

    @QtCore.Slot()
    def insert_new_dpi(self):
        self.dpiInput.setStyleSheet("color: black")

    @QtCore.Slot()
    def update_dpi(self):
        input_text = self.dpiInput.text()
        if input_text == '':
            self.dpiInput.setStyleSheet("color: gray")
            self.dpiInput.setText("300")


# class to handle relaxation time plot
class PlotTau(FigureCanvasQTAgg):
    def __init__(self, parent=None, width=15, height=15, dpi=55):
        super(PlotTau, self).__init__(Figure())

        self.setParent(parent)
        with plt.rc_context({'lines.linewidth': 1, 'lines.linestyle': ':',
                                'xtick.labelsize': 10, 'xtick.major.size': 2,
                                'ytick.labelsize': 10, 'ytick.major.size': 2,
                                'axes.labelsize': 10, 'font.size': 10}):

            self.figure, self.ax = plt.subplots(1, 1, figsize=(width, height), dpi=dpi)
            self.figure.subplots_adjust(left=0.2,bottom=0.2)
            self.canvas = FigureCanvasQTAgg(self.figure)
            FigureCanvasQTAgg.setSizePolicy(self,QtWidgets.QSizePolicy.Expanding,QtWidgets.QSizePolicy.Expanding)
            FigureCanvasQTAgg.updateGeometry(self)
            self.colorbar = None

    def plot(self, mu, tau, T):
        tau = np.transpose(tau) # python - julia compatibility
        sm = plt.cm.ScalarMappable(cmap=cm.viridis, norm=plt.Normalize(vmin=T.min(), vmax=T.max()))
        # single line
        if mu.size == 1:
            self.ax.plot(T,tau,color=tuple(item / 255 for item in gui_color))
            self.ax.set_xlabel(r"$T\ [K]$", fontsize=12)
            self.ax.set_ylabel(r"$\tau$", fontsize=12)
            self.ax.grid(linewidth=0.5)
        # one line for each T
        else:
            sm.set_array([])
            for i in range(T.size):
                self.ax.plot(mu,tau[i,:], color = sm.to_rgba(T[i]))
            if self.colorbar is None:
                self.colorbar = plt.colorbar(sm)
            else:
                self.colorbar.update_normal(sm)
            self.ax.set_xlabel(r"$\mu\ [eV]$", fontsize=12)
            self.ax.set_ylabel(r"$\tau$", fontsize=12)
            self.ax.grid(linewidth=0.5)
        self.draw()


# class to handle transport coefficients plots
class PlotsCanvas(FigureCanvasQTAgg):
    def __init__(self, parent=None, width=5, height=5, dpi=100):
        super(PlotsCanvas, self).__init__(Figure())

        self.setParent(parent)
        self.figure, self.ax = plt.subplots(2, 2, figsize=(width, height), dpi=dpi)
        plt.subplots_adjust(wspace=0.5, hspace=0.5, bottom=0.15)
        self.ax1 = self.ax[0, 0]; self.ax2 = self.ax[0, 1]; self.ax3 = self.ax[1, 0]; self.ax4 = self.ax[1, 1]
        self.canvas = FigureCanvasQTAgg(self.figure)
        FigureCanvasQTAgg.setSizePolicy(self,
                                   QtWidgets.QSizePolicy.Expanding,
                                   QtWidgets.QSizePolicy.Expanding)
        FigureCanvasQTAgg.updateGeometry(self)
        self.point1 = None; self.point2 = None; self.point3 = None; self.point4 = None
        self.colorbar1 = None; self.colorbar2 = None; self.colorbar3 = None; self.colorbar4 = None


    def plot(self, tensor_name, x, y, z, tau_model, exp_data):
        self.figure.suptitle(r"$\tau\ $"+tau_model, y=0.97)
        # plot1: sigma
        if tensor_name == "conductivity":
            # single curve
            if z.size == 1:
                self.ax1.plot(x, np.squeeze(y), "-.", marker='.', fillstyle='none', color='#1f77b4', label="mu=" + str(np.round(z, 4)), zorder=0)
            # one curve for each value of Fermi level
            else:
                sm = plt.cm.ScalarMappable(cmap=cm.viridis, norm=plt.Normalize(vmin=z.min(), vmax=z.max()))
                for i in range(z.size):
                    self.ax1.plot(x, y[i,:], color=sm.to_rgba(z[i]), linewidth=0.5)
                if self.colorbar1 is None:
                    self.colorbar1 = plt.colorbar(sm,ax=self.ax1)
                else:
                    self.colorbar1.update_normal(sm)
            # plot experimental data (if imported)
            if exp_data.df is not None:
                self.ax1.plot(exp_data['temperature'], exp_data[tensor_name], "--", color='dimgray', label="exp")
            self.ax1.ticklabel_format(style="sci", axis='y', scilimits=(3,0))
            # self.ax1.set_xlabel(r"$T\ [K]$")
            self.ax1.set_ylabel(r"$\sigma\ [(\Omega m)^{-1}]$")
            self.ax1.grid(linewidth=0.3)

        # plot2: seebeck
        if tensor_name == "seebeck":
            # single curve
            if z.size == 1:
                self.ax2.plot(x, np.multiply(np.squeeze(y), 1e6), "-.", marker='.', fillstyle='none', color="orange", label="mu=" + str(np.round(z, 4)), zorder=0)
            # one curve for each value of Fermi level
            else:
                sm = plt.cm.ScalarMappable(cmap=cm.viridis, norm=plt.Normalize(vmin=z.min(), vmax=z.max()))
                for i in range(z.size):
                    self.ax2.plot(x, np.multiply(y[i,:], 1e6), color=sm.to_rgba(z[i]), linewidth=0.5)
                if self.colorbar2 is None:
                    self.colorbar2 = plt.colorbar(sm,ax=self.ax2)
                else:
                    self.colorbar2.update_normal(sm)
            if exp_data.df is not None:
                self.ax2.plot(exp_data['temperature'], np.multiply(exp_data[tensor_name], 1e6), "--", color='dimgray', label="exp")
            # self.ax2.set_xlabel(r"$T\ [K]$")
            self.ax2.set_ylabel(r"$S\ [\mu VK^{-1}]$")
            self.ax2.grid(linewidth=0.3)

        # plot3: thermal
        if tensor_name == "thermal":
            # single curve
            if z.size == 1:
                self.ax3.plot(x, np.squeeze(y), "-.", marker='.', fillstyle='none', color="red", label="mu=" + str(np.round(z, 4)), zorder=0)
            else:
                sm = plt.cm.ScalarMappable(cmap=cm.viridis, norm=plt.Normalize(vmin=z.min(), vmax=z.max()))
                for i in range(z.size):
                    self.ax3.plot(x, y[i,:], color=sm.to_rgba(z[i]), linewidth=0.5)
                if self.colorbar3 is None:
                    self.colorbar3 = plt.colorbar(sm,ax=self.ax3)
                else:
                    self.colorbar3.update_normal(sm)
            if exp_data.df is not None:
                self.ax3.plot(exp_data['temperature'], exp_data[tensor_name], "--", color='dimgray', label="exp")
            self.ax3.set_xlabel(r"$T\ [K]$")
            self.ax3.set_ylabel(r"$\kappa_{e}\ [WK^{-1}]$")
            self.ax3.grid(linewidth=0.3)

        # plot4: concentration
        if tensor_name == "concentration":
            # single curve
            if z.size == 1:
                self.ax4.plot(x, np.squeeze(y), "-.", marker='.', fillstyle='none', color="limegreen", label="mu=" + str(np.round(z, 4)), zorder=0)
            # one curve for each value of Fermi level
            else:
                sm = plt.cm.ScalarMappable(cmap=cm.viridis, norm=plt.Normalize(vmin=z.min(), vmax=z.max()))
                for i in range(z.size):
                    self.ax4.plot(x, y[i,:], color=sm.to_rgba(z[i]), linewidth=0.5)
                if self.colorbar4 is None:
                    self.colorbar4 = plt.colorbar(sm,ax=self.ax4)
                else:
                    self.colorbar4.update_normal(sm)
            if exp_data.df is not None:
                self.ax4.plot(exp_data['temperature'], exp_data[tensor_name], "--", color='dimgray', label="exp")
            self.ax4.set_xlabel(r"$T\ [K]$")
            self.ax4.set_ylabel("n")
            self.ax4.grid(linewidth=0.3)

        self.draw()

    # save the plots
    def save(self, fullpath, dpi):
        self.ax1.set_title(""); self.ax2.set_title(""); self.ax3.set_title(""); self.ax4.set_title("")
        plt.savefig(fullpath, dpi=dpi)


if __name__ == "__main__":
    import subprocess

    import platform
    if platform.system() == "Windows":
        import ctypes
        myappid = 'Mstar2t.bonal1l@cmich.edu' # arbitrary string
        ctypes.windll.shell32.SetCurrentProcessExplicitAppUserModelID(myappid)

    # display
    os.environ["QT_AUTO_SCREEN_SCALE_FACTOR"] = "1"
    QtWidgets.QApplication.setAttribute(QtCore.Qt.AA_EnableHighDpiScaling, True)

    # server
    server = subprocess.Popen(['julia', '../run_server.jl'])

    # logo
    app = QtWidgets.QApplication(sys.argv)
    screen = app.primaryScreen()
    loading = LoadingScreen(screen.size())
    app.exec_()

    # GUI
    app = QtWidgets.QApplication.instance()
    app.setStyleSheet(mySetStyleSheet)  # design
    InputWindow = MainWindow(server)
    ui_in = UiInputWindow()
    ui_in.setupUi(InputWindow)
    InputWindow.show()
    ui_in.check_server_status()
    sys.exit(app.exec_())