import sys, os, time, math, datetime, copy, re,  h5py
from collections import OrderedDict
import glob
try:
    from PyQt4.QtCore import *
    from PyQt4.QtGui import *
except:
    from PyQt5.QtCore import *
    from PyQt5.QtGui import *
    from PyQt5.QtWidgets import *
import pyqtgraph as pg
import numpy as np
sys.path.append(os.path.abspath(os.path.realpath(__file__)+'/../../../'))
from SimulationFramework.Modules.read_beam_file import beam as rbfBeam
import SimulationFramework.Modules.read_twiss_file as rtf
sys.path.append(os.path.realpath(__file__)+'/../../../../')

class mainWindow(QMainWindow):
    def __init__(self):
        super(mainWindow, self).__init__()
        self.resize(1800,900)
        self.centralWidget = QWidget()
        self.layout = QVBoxLayout()
        self.centralWidget.setLayout(self.layout)

        self.tab = QTabWidget()
        self.beamPlot = beamPlotWidget()

        self.layout.addWidget(self.beamPlot)

        self.setCentralWidget(self.centralWidget)

        self.setWindowTitle("ASTRA Data Plotter")
        menubar = self.menuBar()
        fileMenu = menubar.addMenu('&File')

        exitAction = QAction('&Exit', self)
        exitAction.setShortcut('Ctrl+Q')
        exitAction.setStatusTip('Exit application')
        exitAction.triggered.connect(self.close)
        fileMenu.addAction(exitAction)

class beamPlotWidget(QWidget):
    # Styles for the plot lines
    colors = [QColor('#F5973A'),QColor('#A95AA1'),QColor('#85C059'),QColor('#0F2080'),QColor('#BDB8AD'), 'r']

    beamParams = OrderedDict([
        ['x', {'quantity': 'x', 'units': 'm', 'name': 'x'}],
        ['y', {'quantity': 'y', 'units': 'm', 'name': 'y'}],
        ['z', {'quantity': 'z', 'units': 'm', 'name': 'z', 'norm': True}],
        ['t', {'quantity': 't', 'units': 's', 'name': 't', 'norm': True}],
        ['cpx', {'quantity': 'cpx', 'units': 'eV', 'name': 'cp<sub>x</sub>'}],
        ['cpy', {'quantity': 'cpy', 'units': 'eV', 'name': 'cp<sub>y</sub>'}],
        ['cpz', {'quantity': 'cpz', 'units': 'eV', 'name': 'cp<sub>z</sub>'}],
    ])

    def __init__(self, **kwargs):
        super(beamPlotWidget, self).__init__(**kwargs)
        ''' These are for reading data files from ASTRA and Elegant '''
        self.beams = {}
        self.twiss = rtf.twiss()

        ''' beamPlotWidget '''
        self.beamPlotWidget = QWidget()
        self.beamPlotLayout = QVBoxLayout()
        self.beamPlotWidget.setLayout(self.beamPlotLayout)
        self.beamPlotPlotWidget = pg.PlotWidget()

        self.beamPlotAxisWidget = QWidget()
        self.beamPlotAxisWidget.setMaximumHeight(100)
        Form = self.beamPlotAxisWidget

        sizePolicy = QSizePolicy(QSizePolicy.Preferred, QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(Form.sizePolicy().hasHeightForWidth())
        Form.setSizePolicy(sizePolicy)
        self.horizontalLayout = QHBoxLayout(Form)
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.beamPlotXAxisCombo = QComboBox(Form)
        self.beamPlotXAxisCombo.addItems(list(self.beamParams.keys()))
        self.beamPlotXAxisCombo.setCurrentIndex(2)
        self.horizontalLayout.addWidget(self.beamPlotXAxisCombo)
        self.beamPlotXAxisNormalise = QCheckBox('Normalise')
        self.beamPlotXAxisNormalise.setChecked(True)
        self.horizontalLayout.addWidget(self.beamPlotXAxisNormalise)
        self.beamPlotYAxisCombo = QComboBox(Form)
        self.beamPlotYAxisCombo.addItems(list(self.beamParams.keys()))
        self.beamPlotYAxisCombo.setCurrentIndex(6)
        self.horizontalLayout.addWidget(self.beamPlotYAxisCombo)
        self.beamPlotYAxisNormalise = QCheckBox('Normalise')
        self.beamPlotYAxisNormalise.setChecked(True)
        self.horizontalLayout.addWidget(self.beamPlotYAxisNormalise)
        self.beamPlotXAxisCombo.currentIndexChanged.connect(self.updateBeamPlot)
        self.beamPlotYAxisCombo.currentIndexChanged.connect(self.updateBeamPlot)
        self.beamPlotXAxisNormalise.stateChanged.connect(self.updateBeamPlot)
        self.beamPlotYAxisNormalise.stateChanged.connect(self.updateBeamPlot)

        self.beamPlotLayout.addWidget(self.beamPlotAxisWidget)
        self.beamPlotLayout.addWidget(self.beamPlotPlotWidget)

        self.layout = QVBoxLayout()
        self.setLayout(self.layout)
        self.layout.addWidget(self.beamPlotWidget)

        ''' used for style cycling '''
        self.plotColor = 0
        self.curves = {}
        self.shadowCurves = []

    def addbeamDataFiles(self, dicts):
        for d in dicts:
            self.addbeamDataFile(**d)

    def addbeamData(self, beamobject, datafile):
        ''' addbeamData - takes a "read_beam_file" object and add a plotItem to the relevant self.curves '''
        ''' add the beam object into the self.beams dictionary '''
        ''' Requires a reference file name '''
        if str(type(beamobject)) == "<class 'SimulationFramework.Modules.read_beam_file.beam'>":
            self.beams[datafile] = beamobject
            color = self.colors[self.plotColor]
            pen = pg.mkBrush(color=color)
            self.curves[datafile] = self.beamPlotPlotWidget.plot([], pen=None,
            symbolBrush=pen, symbolSize=5, symbolPen=None)
            self.curves[datafile].sigClicked.connect(lambda: self.highlightPlot(datafile))
            self.plotColor += 1
        self.updateBeamPlot()

    def addbeamDataFile(self, directory, filename):
        ''' addbeamDataFile - read the data file and add a plotItem to the relevant self.curves '''
        ''' load the data file into the self.beams dictionary '''
        datafile = directory + '/' + filename
        if os.path.isfile(datafile):
            beam = rbfBeam()
            beam.read_HDF5_beam_file(datafile)
            self.addbeamData(beam, datafile)

    def updateBeamPlot(self):
        xdict = self.beamParams[str(self.beamPlotXAxisCombo.currentText())]
        ydict = self.beamParams[str(self.beamPlotYAxisCombo.currentText())]
        for d in self.curves:
            x = getattr(self.beams[d], str(self.beamPlotXAxisCombo.currentText()))
            if self.beamPlotXAxisNormalise.isChecked():
                x = x - np.mean(x)
            y = getattr(self.beams[d], str(self.beamPlotYAxisCombo.currentText()))
            if self.beamPlotYAxisNormalise.isChecked():
                y = y - np.mean(y)
            self.curves[d].setData(x=x, y=y)
        self.beamPlotPlotWidget.setLabel('left', text=ydict['name'], units=ydict['units'])
        self.beamPlotPlotWidget.setLabel('bottom', text=xdict['name'], units=xdict['units'])

    def removePlot(self, directory):
        ''' finds all beam plots based on a directory name, and removes them '''
        if not isinstance(directory, (list, tuple)):
            directory = [directory]
        indexname = directory[-1]
        if directory in self.beamPlotItems:
            for entry in self.beamplotLayout:
                if entry == 'next_row':
                    pass
                else:
                    self.beamPlotWidgets[entry['label']].removeItem(self.beamPlotItems[indexname][entry['label']])

    def highlightPlot(self, name):
        ''' highlights a particular plot '''
        # print('highligher clicked! = ', name)
        if not isinstance(name, (list, tuple)):
            name = [name]
        for n in name:
            self.addShadowPen(n)
        for n in self.curves.keys():
            if n in self.shadowCurves or not len(self.shadowCurves) > 0:
                self.setPenAlpha(n, 255, 3)
            else:
                self.setPenAlpha(n, 10, 3)

    def addShadowPen(self, name):
        # curve = self.curves[name]
        if not name in self.shadowCurves:
            self.shadowCurves.append(name)
            # pen = curve.opts['symbolBrush']
            # shadowpencolor = pen.color()
            # shadowpencolor.setAlpha(100)
            # shadowpen = pg.mkPen(color=shadowpencolor, width=6)
            # curve.setSymbolPen(shadowpen)
        else:
            self.shadowCurves.remove(name)
            # curve.setSymbolPen(None)
            # curve.opts['symbolPen'] = None

    def setPenAlpha(self, name, alpha=255, width=3):
        curve = self.curves[name]
        pen = curve.opts['symbolBrush']
        pencolor = pen.color()
        pencolor.setAlpha(alpha)
        pen = pg.mkBrush(color=pencolor, width=width, style=pen.style())
        curve.setSymbolBrush(pen)

def main():
    app = QApplication(sys.argv)
    # pg.setConfigOptions(antialias=True)
    pg.setConfigOption('background', 'w')
    pg.setConfigOption('foreground', 'k')
    ex = mainWindow()
    ex.show()
    ex.beamPlot.addbeamDataFiles([
    {'directory': 'OnlineModel_test_data/basefiles_4_250pC', 'filename': 'CLA-S02-APER-01.hdf5'},
    {'directory': 'OnlineModel_test_data/test_4', 'filename': 'CLA-S07-APER-01.hdf5'}])
    sys.exit(app.exec_())

if __name__ == '__main__':
   main()
