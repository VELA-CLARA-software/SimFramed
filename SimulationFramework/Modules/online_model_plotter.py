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
from pyqtgraph import LegendItem, mkPen, mkBrush, LabelItem, TableWidget, GraphicsLayoutWidget, setConfigOption, \
setConfigOptions, InfiniteLine, ImageItem, GraphicsView, GraphicsLayout, AxisItem, ViewBox, PlotDataItem, colorStr, mkColor, ImageView, PlotItem
from pyqtgraph.graphicsItems.LegendItem import ItemSample
import argparse
import numpy as np
sys.path.append(os.path.abspath(os.path.realpath(__file__)+'/../../../'))
# print (sys.path)
import SimulationFramework.Modules.read_beam_file as raf
import SimulationFramework.Modules.read_twiss_file as rtf
from SimulationFramework.Modules.online_model_twissPlot import twissPlotWidget
from SimulationFramework.Modules.online_model_slicePlot import slicePlotWidget
sys.path.append(os.path.realpath(__file__)+'/../../../../')

class mainWindow(QMainWindow):
    def __init__(self, parent = None, directory='.'):
        super(mainWindow, self).__init__(parent)
        self.resize(1800,900)

        self.setWindowTitle("ASTRA Data Plotter")
        menubar = self.menuBar()
        fileMenu = menubar.addMenu('&File')

        # reloadSettingsAction = QAction('Reload Settings', self)
        # reloadSettingsAction.setStatusTip('Reload Settings YAML File')
        # reloadSettingsAction.triggered.connect(self.picklePlot.reloadSettings)
        # fileMenu.addAction(reloadSettingsAction)

        exitAction = QAction('&Exit', self)
        exitAction.setShortcut('Ctrl+Q')
        exitAction.setStatusTip('Exit application')
        exitAction.triggered.connect(self.close)
        fileMenu.addAction(exitAction)


        self.onlineModelPlotter = onlineModelPlotter(directory=directory)
        self.setCentralWidget(self.onlineModelPlotter)


class onlineModelPlotter(QWidget):
    def __init__(self, parent = None, directory='.'):
        super(onlineModelPlotter, self).__init__()
        self.directory = directory

        self.layout = QVBoxLayout()
        self.setLayout(self.layout)

        self.tab = QTabWidget()

        self.tabWidget = QTabWidget()

        self.twissPlotWidget = twissPlotWidget()
        self.slicePlotWidget = slicePlotWidget()

        self.folderButton = QPushButton('Select Directory')
        self.folderLineEdit = QLineEdit()
        # self.folderLineEdit.setReadOnly(True)
        self.folderLineEdit.setText(self.directory)
        self.reloadButton = QPushButton()
        self.reloadButton.setIcon(qApp.style().standardIcon(QStyle.SP_BrowserReload	))
        self.folderWidget = QGroupBox()
        self.folderLayout = QHBoxLayout()
        self.folderLayout.addWidget(self.folderButton)
        self.folderLayout.addWidget(self.folderLineEdit)
        self.folderLayout.addWidget(self.reloadButton)
        self.folderWidget.setLayout(self.folderLayout)
        self.folderWidget.setMaximumWidth(800)
        self.reloadButton.clicked.connect(lambda: self.changeDirectory(self.directory))
        self.folderButton.clicked.connect(self.changeDirectory)

        self.fileSelector = QComboBox()
        self.fileSelector.currentIndexChanged.connect(self.updateScreenCombo)
        self.screenSelector = QComboBox()
        self.screenSelector.currentIndexChanged.connect(self.updatePlotButtons)
        self.plotScreenButton = QPushButton('Plot')
        self.plotScreenButton.clicked.connect(self.loadScreen)
        self.removeScreenButton = QPushButton('Remove')
        self.removeScreenButton.clicked.connect(self.removeScreen)
        self.removeScreenButton.setEnabled(False)
        self.beamWidget = QGroupBox()
        self.beamLayout = QHBoxLayout()
        self.beamLayout.addWidget(self.fileSelector)
        self.beamLayout.addWidget(self.screenSelector)
        self.beamLayout.addWidget(self.plotScreenButton)
        self.beamLayout.addWidget(self.removeScreenButton)
        self.beamWidget.setLayout(self.beamLayout)
        self.beamWidget.setMaximumWidth(800)
        self.beamWidget.setVisible(False)

        self.sliceWidget = QGroupBox()
        self.sliceLayout = QHBoxLayout()
        self.sliceWidget.setLayout(self.sliceLayout)
        self.sliceLayout.addWidget(self.slicePlotWidget.slicePlotSliceWidthWidget)
        self.sliceWidget.setVisible(False)
        self.sliceWidget.setMaximumWidth(150)

        self.folderBeamWidget = QWidget()
        self.folderBeamLayout = QHBoxLayout()
        self.folderBeamLayout.setAlignment(Qt.AlignLeft);
        self.folderBeamWidget.setLayout(self.folderBeamLayout)
        self.folderBeamLayout.addWidget(self.folderWidget)
        self.folderBeamLayout.addWidget(self.beamWidget)
        self.folderBeamLayout.addWidget(self.sliceWidget)

        self.tabWidget.addTab(self.twissPlotWidget,'Twiss Plots')
        # self.tabWidget.addTab(self.beamPlotWidget,'Beam Plots')
        self.tabWidget.addTab(self.slicePlotWidget,'Slice Beam Plots')
        self.tabWidget.currentChanged.connect(self.changeTab)
        self.layout.addWidget(self.folderBeamWidget)
        self.layout.addWidget(self.tabWidget)

        self.plotType = 'Twiss'
        self.twissDataCounter = 0
        self.changeDirectory(self.directory)

    def loadTwissDataFile(self):
        # if self.plotType == 'Twiss':
        self.twissPlotWidget.addTwissDirectory([{'directory': self.directory, 'sections': 'All'}], name=self.twissDataCounter)
        # elif self.plotType == 'Beam' or self.plotType == 'Slice':
        # if hasattr(self,'beamFileName') and os.path.isfile(self.directory+'/'+self.beamFileName):
        #     # starttime = time.time()
        #     self.beam.read_HDF5_beam_file(self.directory+'/'+self.beamFileName)
        #     # print 'reading file took ', time.time()-starttime, 's'
        #     # print 'Read file: ', self.beamFileName
        #     if self.plotType == 'Beam':
        #         self.plotDataBeam()
        #     else:
        #         self.changeSliceLength()
        #         # self.plotDataSlice()
        self.twissDataCounter += 1

    def loadBeamDataFile(self):
        if hasattr(self,'beamFileName') and os.path.isfile(self.directory+'/'+self.beamFileName):
            self.slicePlotWidget.addsliceDataFile(self.directory+'/'+self.beamFileName)
        self.updatePlotButtons()

    def updatePlotButtons(self):
        if len(self.screenpositions) > 0 and not self.fileSelector.currentText() == '':
            run = self.screenpositions[str(self.fileSelector.currentText())]['run']
            self.beamFileName = str(self.screenpositions[str(self.fileSelector.currentText())]['screennames'][self.screenSelector.currentIndex()])+'.hdf5'
            if self.directory+'/'+self.beamFileName in self.slicePlotWidget.curves:
                self.plotScreenButton.setText('Update')
                self.removeScreenButton.setEnabled(True)
            else:
                self.plotScreenButton.setText('Plot')
                self.removeScreenButton.setEnabled(False)

    def changeTab(self, i):
        if self.tabWidget.tabText(i) == 'Beam Plots':
            self.plotType = 'Beam'
            self.beamWidget.setVisible(True)
            self.sliceWidget.setVisible(True)
        elif self.tabWidget.tabText(i) == 'Slice Beam Plots':
            self.plotType = 'Slice'
            self.beamWidget.setVisible(True)
            self.sliceWidget.setVisible(True)
        else:
            self.plotType = 'Twiss'
            self.beamWidget.setVisible(False)
            self.sliceWidget.setVisible(False)
        # self.loadDataFile()

    def changeDirectory(self, directory=None):
        try:
            self.fileSelector.currentIndexChanged.disconnect(self.updateScreenCombo)
        except:
            pass
        if directory == None or directory == False:
            self.directory = str(QFileDialog.getExistingDirectory(self, "Select Directory", self.directory, QFileDialog.ShowDirsOnly))
        else:
            self.directory = directory
        self.folderLineEdit.setText(self.directory)
        self.currentFileText = self.fileSelector.currentText()
        self.currentScreenText = self.screenSelector.currentText()
        self.getScreenFiles()
        self.updateFileCombo()
        self.updateScreenCombo()
        self.loadTwissDataFile()
        self.fileSelector.currentIndexChanged.connect(self.updateScreenCombo)

    def getSposition(self, file):
        file = h5py.File(self.directory+'/'+file+'.hdf5', "r")
        zpos = file.get('/Parameters/Starting_Position')[2]
        if abs(zpos) < 0.01:
            print(zpos, file)
        return zpos

    def getScreenFiles(self):
        self.screenpositions = {}
        files = glob.glob(self.directory+'/*.hdf5')
        filenames = ['-'.join(os.path.basename(f).split('-')[:2]) for f in files]
        runnumber = ['001' for f in filenames]
        for f in filenames:
            files = glob.glob(self.directory+'/'+f+'*.hdf5')
            screennames = sorted([os.path.basename(s).split('.')[0] for s in files], key=lambda x: self.getSposition(x))
            screenpositions = [self.getSposition(s) for s in screennames]
            self.screenpositions[f] = {'screenpositions': screenpositions, 'screennames': screennames, 'run': '001'}

    def updateFileCombo(self):
        self.fileSelector.clear()
        i = -1
        screenfirstpos = []
        for f in self.screenpositions:
            if len(self.screenpositions[f]['screenpositions']) > 0:
                screenfirstpos.append([f, min(self.screenpositions[f]['screenpositions'])])
        screenfirstpos = np.array(screenfirstpos)
        sortedscreennames = sorted(screenfirstpos, key=lambda x: float(x[1]))
        for f in sortedscreennames:
            self.fileSelector.addItem(f[0])
            i += 1
            if f[0] == self.currentFileText:
                self.fileSelector.setCurrentIndex(i)

    def loadScreen(self, i):
        if len(self.screenpositions) > 0 and not self.fileSelector.currentText() == '':
            run = self.screenpositions[str(self.fileSelector.currentText())]['run']
            self.beamFileName = str(self.screenpositions[str(self.fileSelector.currentText())]['screennames'][self.screenSelector.currentIndex()])+'.hdf5'
            self.loadBeamDataFile()

    def removeScreen(self, i):
        if len(self.screenpositions) > 0 and not self.fileSelector.currentText() == '':
            self.slicePlotWidget.removeCurve(self.directory+'/'+self.beamFileName)
        self.updatePlotButtons()

    def updateScreenCombo(self):
        self.screenSelector.clear()
        i = -1
        if len(self.screenpositions) > 0 and not self.fileSelector.currentText() == '':
            for i, s in enumerate(self.screenpositions[str(self.fileSelector.currentText())]['screennames']):
                self.screenSelector.addItem(str(s)+' ('+str(self.screenpositions[str(self.fileSelector.currentText())]['screenpositions'][i])+')')
                i += 1
                if s == self.currentScreenText:
                    self.screenSelector.setCurrentIndex(i)

setConfigOptions(antialias=True)
setConfigOption('background', 'w')
setConfigOption('foreground', 'k')
def main():
    # global app
    # args = parser.parse_args()
    app = QApplication(sys.argv)
    setConfigOptions(antialias=True)
    setConfigOption('background', 'w')
    setConfigOption('foreground', 'k')
    ex = mainWindow(directory='.')
    ex.show()
    sys.exit(app.exec_())

if __name__ == '__main__':
   main()
