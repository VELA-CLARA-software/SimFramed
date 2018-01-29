import sys, os, time, math, datetime, copy
import glob
from PyQt4.QtCore import QObject, pyqtSignal, QThread, QTimer, QRectF
from PyQt4.QtGui import QApplication, QMainWindow, QWidget, QPushButton, QVBoxLayout, QHBoxLayout, QComboBox, QTabWidget, QLineEdit, QFileDialog, QLabel, QAction, QPixmap
from pyqtgraph import LegendItem, mkPen, mkBrush, LabelItem, TableWidget, GraphicsLayoutWidget, setConfigOption, setConfigOptions, InfiniteLine, ImageItem
from pyqtgraph.graphicsItems.LegendItem import ItemSample
import argparse
import imageio
import numpy as np
import read_beam_file as raf
import read_twiss_file as rtf

beam = raf.beam()
twiss = rtf.twiss()

parser = argparse.ArgumentParser(description='Plot ASTRA Data Files')
parser.add_argument('-d', '--directory', default='.')


class mainWindow(QMainWindow):
    def __init__(self, parent = None, directory='.'):
        super(mainWindow, self).__init__(parent)
        global app
        self.resize(1800,900)
        self.centralWidget = QWidget()
        self.layout = QHBoxLayout()
        self.centralWidget.setLayout(self.layout)

        self.tab = QTabWidget()
        self.astraPlot = astraPlotWidget(directory)

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

class astraPlotWidget(QWidget):

    twissLayout = [{'name': 'sigma_x', 'range': [0,1], 'scale': 1e3},
                   {'name': 'sigma_y', 'range': [0,1], 'scale': 1e3},
                   {'name': 'kinetic_energy', 'range': [0,250], 'scale': 1e-6},
                   'next_row',
                   {'name': 'sigma_p', 'range': [0,0.015], 'scale': 1e6},
                   {'name': 'sigma_z', 'range': [0,0.6], 'scale': 1e3},
                   {'name': 'enx', 'range': [0.5,1.5], 'scale': 1e6},
                   'next_row',
                   {'name': 'eny', 'range':  [0.5,1.5], 'scale': 1e6},
                   {'name': 'beta_x', 'range': [0,150], 'scale': 1},
                   {'name': 'beta_y', 'range': [0,150], 'scale': 1},
                  ]

    def __init__(self, directory='.', **kwargs):
        super(astraPlotWidget, self).__init__(**kwargs)

        self.plotWidget = GraphicsLayoutWidget()
        self.layout = QVBoxLayout()
        self.setLayout(self.layout)
        self.layout.addWidget(self.plotWidget)
        self.twissPlot = True
        self.latticePlotData = imageio.imread('../lattice_plot.png')
        self.loadDataFile(directory)

    def loadDataFile(self, directory):
        if self.twissPlot:
            files = sorted(glob.glob(directory+"/*Xemit*"))
            twiss.read_astra_emit_files(files)
            self.plotDataTwiss()
        else:
            beam.read_astra_beam_file(filename)
            self.twissFile = False
            self.beamFile = True

    def plotDataTwiss(self):
        for entry in self.twissLayout:
            if entry == 'next_row':
                self.plotWidget.nextRow()
            else:
                p = self.plotWidget.addPlot(title=entry['name'])
                p.showGrid(x=True, y=True)
                p.vb.setYRange(*entry['range'])
                latticePlot = ImageItem(self.latticePlotData)
                latticePlot.setOpts(axisOrder='row-major')
                p.vb.addItem(latticePlot)
                latticePlot.setZValue(-1)  # make sure this image is on top
                latticePlot.setOpacity(0.5)
                scaleY = 0.05*(entry['range'][1] - entry['range'][0])
                rect = QRectF(0, entry['range'][0] + scaleY, twiss['z'][-1], 4*scaleY)
                latticePlot.setRect(rect)
                # latticePlot.scale(scaleX, scaleY)
                x = twiss['z']
                y = twiss[entry['name']]*entry['scale']
                plot = p.plot(x=x, y=y, pen=mkPen('b', width=3))

def main():
    global app
    args = parser.parse_args()
    app = QApplication(sys.argv)
    setConfigOptions(antialias=True)
    setConfigOption('background', 'w')
    setConfigOption('foreground', 'k')
    # app.setStyle(QStyleFactory.create("plastique"))
    ex = mainWindow(directory=args.directory)
    ex.show()
    sys.exit(app.exec_())

if __name__ == '__main__':
   main()
