import sys, os, time, math, datetime, copy, re
import glob
from PyQt4.QtCore import QObject, pyqtSignal, QThread, QTimer, QRectF, Qt
from PyQt4.QtGui import QApplication, QMainWindow, QWidget, QPushButton, QVBoxLayout, QHBoxLayout, QComboBox, QTabWidget, QLineEdit, QFileDialog, QLabel, QAction, QPixmap, qApp, QStyle, QGroupBox
from pyqtgraph import LegendItem, mkPen, mkBrush, LabelItem, TableWidget, GraphicsLayoutWidget, setConfigOption, setConfigOptions, InfiniteLine, ImageItem, GraphicsView, GraphicsLayout
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
        self.directory = directory
        self.resize(1800,900)
        self.centralWidget = QWidget()
        self.layout = QVBoxLayout()
        self.centralWidget.setLayout(self.layout)

        self.tab = QTabWidget()
        self.astraPlot = astraPlotWidget(self.directory)

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

    twissplotLayout = [{'name': 'sigma_x', 'range': [0,1], 'scale': 1e3},
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
        self.directory = directory
        self.twissPlotView = GraphicsView(useOpenGL=True)
        self.twissPlotWidget = GraphicsLayout()
        self.twissPlotView.setCentralItem(self.twissPlotWidget)

        self.latticePlotData = imageio.imread('lattice_plot.png')
        self.latticePlots = {}
        self.twissPlots = {}
        i = -1
        for entry in self.twissplotLayout:
            if entry == 'next_row':
                self.twissPlotWidget.nextRow()
            else:
                i += 1
                p = self.twissPlotWidget.addPlot(title=entry['name'])
                p.showGrid(x=True, y=True)
                vb = p.vb
                vb.setYRange(*entry['range'])
                latticePlot = ImageItem(self.latticePlotData)
                latticePlot.setOpts(axisOrder='row-major')
                vb.addItem(latticePlot)
                latticePlot.setZValue(-1)  # make sure this image is on top
                latticePlot.setOpacity(0.5)
                self.twissPlots[entry['name']] = p.plot(pen=mkPen('b', width=3))
                self.latticePlots[p.vb] = latticePlot
                p.vb.sigRangeChanged.connect(self.scaleLattice)
        self.beamPlotWidget = QWidget()
        self.beamPlotLayout = QVBoxLayout()
        self.beamPlotWidget.setLayout(self.beamPlotLayout)
        self.beamPlotView = GraphicsView(useOpenGL=True)
        self.beamPlotWidgetGraphicsLayout = GraphicsLayout()
        p = self.beamPlotWidgetGraphicsLayout.addPlot(title='beam')
        p.showGrid(x=True, y=True)
        self.beamPlot = p.plot(pen=None, symbol='+')
        self.beamPlotView.setCentralItem(self.beamPlotWidgetGraphicsLayout)
        self.beamPlotXAxisCombo = QComboBox()
        self.beamPlotXAxisCombo.addItems(['x','y','zn','cpx','cpy','BetaGamma'])
        self.beamPlotYAxisCombo = QComboBox()
        self.beamPlotYAxisCombo.addItems(['x','y','zn','cpx','cpy','BetaGamma'])
        self.beamPlotAxisWidget = QWidget()
        self.beamPlotAxisLayout = QHBoxLayout()
        self.beamPlotAxisWidget.setLayout(self.beamPlotAxisLayout)
        self.beamPlotAxisLayout.addWidget(self.beamPlotXAxisCombo)
        self.beamPlotAxisLayout.addWidget(self.beamPlotYAxisCombo)
        self.beamPlotXAxisCombo.currentIndexChanged.connect(self.plotDataBeam)
        self.beamPlotYAxisCombo.currentIndexChanged.connect(self.plotDataBeam)
        # self.beamPlotXAxisCombo.setCurrentIndex(2)
        # self.beamPlotYAxisCombo.setCurrentIndex(5)
        self.beamPlotLayout.addWidget(self.beamPlotAxisWidget)
        self.beamPlotLayout.addWidget(self.beamPlotView)

        self.layout = QVBoxLayout()
        self.setLayout(self.layout)

        self.tabWidget = QTabWidget()

        self.folderButton = QPushButton('Select Directory')
        self.folderLineEdit = QLineEdit()
        self.folderLineEdit.setReadOnly(True)
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
        self.reloadButton.clicked.connect(self.loadDataFile)
        self.folderButton.clicked.connect(self.changeDirectory)

        self.fileSelector = QComboBox()
        self.fileSelector.currentIndexChanged.connect(self.updateScreenCombo)
        self.screenSelector = QComboBox()
        self.screenSelector.currentIndexChanged.connect(self.changeScreen)
        self.beamWidget = QGroupBox()
        self.beamLayout = QHBoxLayout()
        self.beamLayout.addWidget(self.fileSelector)
        self.beamLayout.addWidget(self.screenSelector)
        self.beamWidget.setLayout(self.beamLayout)
        self.beamWidget.setMaximumWidth(800)
        self.beamWidget.setVisible(False)

        self.folderBeamWidget = QWidget()
        self.folderBeamLayout = QHBoxLayout()
        self.folderBeamLayout.setAlignment(Qt.AlignLeft);
        self.folderBeamWidget.setLayout(self.folderBeamLayout)
        self.folderBeamLayout.addWidget(self.folderWidget)
        self.folderBeamLayout.addWidget(self.beamWidget)

        self.tabWidget.addTab(self.twissPlotView,'Twiss Plots')
        self.tabWidget.addTab(self.beamPlotWidget,'Beam Plots')
        self.tabWidget.currentChanged.connect(self.changeTab)
        self.layout.addWidget(self.folderBeamWidget)
        self.layout.addWidget(self.tabWidget)

        self.twissPlot = True
        self.changeDirectory(self.directory)

    def changeTab(self, i):
        if self.tabWidget.tabText(i) == 'Beam Plots':
            self.twissPlot = False
            self.beamWidget.setVisible(True)
        else:
            self.twissPlot = True
            self.beamWidget.setVisible(False)
        self.loadDataFile()

    def changeDirectory(self, directory=None):
        if directory == None:
            self.directory = str(QFileDialog.getExistingDirectory(self, "Select Directory", self.directory, QFileDialog.ShowDirsOnly))
        else:
            self.directory = directory
        self.folderLineEdit.setText(self.directory)
        self.currentFileText = self.fileSelector.currentText()
        self.currentScreenText = self.screenSelector.currentText()
        self.getScreenFiles()
        self.updateFileCombo()
        self.updateScreenCombo()
        self.loadDataFile()

    def getScreenFiles(self):
        self.screenpositions = {}
        files = glob.glob(self.directory+'/*.????.???')
        filenames = ['.'.join(os.path.basename(f).split('.')[:2]) for f in files]
        runnumber = [os.path.basename(f).split('.')[-1] for f in files]
        for f, r in list(set(zip(filenames, runnumber))):
            files = glob.glob(self.directory+'/'+f+'.????.???')
            screenpositions = [re.search('\d\d\d\d', s).group(0) for s in files]
            self.screenpositions[f] = {'screenpositions': screenpositions, 'run': r}

    def updateFileCombo(self):
        self.fileSelector.clear()
        i = -1
        for f in self.screenpositions:
            self.fileSelector.addItem(f)
            i += 1
            if f == self.currentFileText:
                self.fileSelector.setCurrentIndex(i)


    def changeScreen(self, i):
        run = self.screenpositions[str(self.fileSelector.currentText())]['run']
        self.beamFileName = str(self.fileSelector.currentText())+'.'+str(self.screenSelector.currentText())+'.'+str(run)
        # print 'beamFileName = ', self.beamFileName
        self.loadDataFile()

    def updateScreenCombo(self):
        self.screenSelector.clear()
        i = -1
        for s in self.screenpositions[str(self.fileSelector.currentText())]['screenpositions']:
            self.screenSelector.addItem(s)
            i += 1
            if s == self.currentScreenText:
                self.screenSelector.setCurrentIndex(i)

    def loadDataFile(self):
        if self.twissPlot:
            files = sorted(glob.glob(self.directory+"/*Xemit*"))
            twiss.read_astra_emit_files(files)
            self.plotDataTwiss()
        else:
            if hasattr(self,'beamFileName') and os.path.isfile(self.directory+'/'+self.beamFileName):
                beam.read_astra_beam_file(self.directory+'/'+self.beamFileName)
                self.plotDataBeam()

    def plotDataTwiss(self):
        # self.latticePlots = {}
        # self.twissPlotWidget.clear()
        # i = -1
        for entry in self.twissplotLayout:
            if entry == 'next_row':
                pass
            else:
                x = twiss['z']
                y = twiss[entry['name']]*entry['scale']
                self.twissPlots[entry['name']].setData(x=x, y=y, pen=mkPen('b', width=3))

    def plotDataBeam(self):
        self.beamPlot.setData(x=getattr(beam, str(self.beamPlotXAxisCombo.currentText())), y=getattr(beam, str(self.beamPlotYAxisCombo.currentText())), pen=None, symbol='+')

    def scaleLattice(self, vb, range):
        yrange = range[1]
        scaleY = 0.05*abs(yrange[1] - yrange[0])
        rect = QRectF(0, yrange[0] + 2*scaleY, twiss['z'][-1], 4*scaleY)
        self.latticePlots[vb].setRect(rect)

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
