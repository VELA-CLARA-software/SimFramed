import sys, os, time, math, datetime, copy, h5py
import argparse
from PyQt4.QtGui import QApplication, QMainWindow, QWidget, QPushButton, QVBoxLayout, QHBoxLayout, QComboBox, QTabWidget, QLineEdit, QFileDialog, QLabel, QAction, QPixmap, qApp, QStyle, QGroupBox
from pyqtgraph import LegendItem, mkPen, mkBrush, LabelItem, TableWidget, GraphicsLayoutWidget, setConfigOption, setConfigOptions, InfiniteLine, ImageItem, GraphicsView, GraphicsLayout

import astra_plotter_test as astraplotter
import numpy as np
import read_beam_file as raf
import read_twiss_file as rtf

parser = argparse.ArgumentParser(description='Plot ASTRA Data Files')
parser.add_argument('filename')


class mainWindow(QMainWindow):
    def __init__(self, app = None, filename='.'):
        super(mainWindow, self).__init__()
        self.resize(1800,900)
        self.centralWidget = QWidget()
        self.layout = QHBoxLayout()
        self.centralWidget.setLayout(self.layout)

        self.tab = QTabWidget()
        self.astraPlot = astraHDFPlotWidget(filename)

        self.layout.addWidget(self.astraPlot)

        self.setCentralWidget(self.centralWidget)

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
        exitAction.triggered.connect(app.quit)
        fileMenu.addAction(exitAction)

class astraHDFPlotWidget(astraplotter.astraPlotWidget):

    def __init__(self, filename=None, **kwargs):
        super(astraHDFPlotWidget, self).__init__(**kwargs)
        print 'filename = ', filename
        self.folderButton.setText('Select File')
        self.filename = filename
        print 'self.filename = ', self.filename
        self.h5file = h5py.File(self.filename, "r")
        self.changeDirectory(filename)

    def changeDirectory(self, filename=None):
        if filename == None or filename == False:
            self.filename = str(QFileDialog.getOpenFileName(self, 'Select HDF file', '.', filter="HDF5 files (*.h5 *.hdf5);;", selectedFilter="HDF5 files (*.h5 *.hdf5)"))
        else:
            self.filename = filename
        self.h5file = h5py.File(self.filename, "r")
        self.folderLineEdit.setText(self.filename)
        self.currentFileText = self.fileSelector.currentText()
        self.currentScreenText = self.screenSelector.currentText()
        self.getScreenFiles()
        self.updateFileCombo()
        self.updateScreenCombo()
        self.loadDataFile()

    def getScreenFiles(self):
        try:
            self.screenpositions = {}
            h5screens = self.h5file.get('screens')
            screenpositions = h5screens.keys()
            self.screenpositions[self.filename] = {'screenpositions': screenpositions}
        except:
            pass

    def changeScreen(self, i):
        try:
            h5screens = self.h5file.get('screens')
            self.beamData = np.array(h5screens.get(str(self.screenSelector.currentText())))
            self.loadDataFile()
        except:
            pass

    def loadDataFile(self):
        if os.path.isfile(self.filename):
            if self.plotType == 'Twiss':
                self.twiss.read_hdf_summary(self.filename)
                self.plotDataTwiss()
            elif self.plotType == 'Beam' or self.plotType == 'Slice':
                self.beam.read_hdf5_beam(self.beamData)
                if self.plotType == 'Beam':
                    self.plotDataBeam()
                else:
                    self.beam.bin_time()
                    self.plotDataSlice()

def main():
    args = parser.parse_args()
    app = QApplication(sys.argv)
    setConfigOptions(antialias=True)
    setConfigOption('background', 'w')
    setConfigOption('foreground', 'k')
    # app.setStyle(QStyleFactory.create("plastique"))
    ex = mainWindow(app, filename=args.filename)
    ex.show()
    sys.exit(app.exec_())

if __name__ == '__main__':
   main()
