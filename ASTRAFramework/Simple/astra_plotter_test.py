import sys, os, time, math, datetime, copy, re
import glob
from PyQt4.QtCore import QObject, pyqtSignal, QThread, QTimer, QRectF, Qt
from PyQt4.QtGui import * #QApplication, QMainWindow, QWidget, QPushButton, QVBoxLayout, QHBoxLayout, QComboBox, QTabWidget, QLineEdit, QFileDialog, QLabel, QAction, QPixmap, qApp, QStyle, QGroupBox, QSpinBox
from pyqtgraph import LegendItem, mkPen, mkBrush, LabelItem, TableWidget, GraphicsLayoutWidget, setConfigOption, \
setConfigOptions, InfiniteLine, ImageItem, GraphicsView, GraphicsLayout, AxisItem, ViewBox, PlotDataItem, colorStr, mkColor
from pyqtgraph.graphicsItems.LegendItem import ItemSample
import argparse
import imageio
import numpy as np
import read_beam_file as raf
import read_twiss_file as rtf

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
        self.beam = raf.beam()
        self.twiss = rtf.twiss()
        self.directory = directory

        ''' twissPlotWidget '''
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

        ''' beamPlotWidget '''
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

        ''' slicePlotWidget '''
        self.sliceParams = [{'name': 'slice_normalized_horizontal_emittance', 'units': 'm-rad', 'text': 'enx'},
        {'name': 'slice_normalized_vertical_emittance', 'units': 'm-rad', 'text': 'eny'},
        {'name': 'slice_peak_current', 'units': 'A', 'text': 'PeakI'},
        {'name': 'slice_relative_momentum_spread', 'units': '%', 'text': 'sigma-p'},
        ]
        self.slicePlotWidget = QWidget()
        self.slicePlotLayout = QVBoxLayout()
        self.slicePlotWidget.setLayout(self.slicePlotLayout)
        # self.slicePlotView = GraphicsView(useOpenGL=True)
        self.slicePlotWidgetGraphicsLayout = GraphicsLayoutWidget()
        # self.slicePlots = {}
        self.slicePlotCheckbox = {}
        self.curve = {}
        self.sliceaxis = {}
        self.slicePlotCheckboxWidget = QWidget()
        self.slicePlotCheckboxLayout = QVBoxLayout()
        self.slicePlotCheckboxWidget.setLayout(self.slicePlotCheckboxLayout)
        self.slicePlot = self.slicePlotWidgetGraphicsLayout.addPlot(title='Slice',row=0,col=50)
        self.slicePlot.showAxis('left', False)
        self.slicePlot.showGrid(x=True, y=True)
        i = -1;
        colors = ['b','r','g','k']
        for param in self.sliceParams:
            i += 1;
            axis = AxisItem("left")
            labelStyle = {'color': '#'+colorStr(mkColor(colors[i]))[0:-2]}
            axis.setLabel(text=param['text'], units=param['units'],**labelStyle)
            viewbox = ViewBox()
            axis.linkToView(viewbox)
            viewbox.setXLink(self.slicePlot.vb)
            self.sliceaxis[param['name']] = [axis, viewbox]
            self.curve[param['name']] = PlotDataItem(pen=colors[i], symbol='+')
            viewbox.addItem(self.curve[param['name']])
            col = self.findFirstEmptyColumnInGraphicsLayout()
            self.slicePlotWidgetGraphicsLayout.ci.addItem(axis, row = 0, col = col,  rowspan=1, colspan=1)
            self.slicePlotWidgetGraphicsLayout.ci.addItem(viewbox, row=0, col=50)
            p.showGrid(x=True, y=True)
            # self.slicePlots[param] = self.slicePlot.plot(pen=colors[i], symbol='+')
            self.slicePlotCheckbox[param['name']] = QCheckBox(param['text'])
            self.slicePlotCheckboxLayout.addWidget(self.slicePlotCheckbox[param['name']])
            self.slicePlotCheckbox[param['name']].stateChanged.connect(self.plotDataSlice)
        # self.slicePlotView.setCentralItem(self.slicePlotWidgetGraphicsLayout)
        self.slicePlotSliceWidthWidget = QSpinBox()
        self.slicePlotSliceWidthWidget.setMaximum(1000)
        self.slicePlotSliceWidthWidget.setValue(100)
        self.slicePlotSliceWidthWidget.setSingleStep(10)
        self.slicePlotSliceWidthWidget.setSuffix("fs")
        self.slicePlotSliceWidthWidget.setSpecialValueText('Automatic')
        self.slicePlotAxisWidget = QWidget()
        self.slicePlotAxisLayout = QHBoxLayout()
        self.slicePlotAxisWidget.setLayout(self.slicePlotAxisLayout)
        self.slicePlotAxisLayout.addWidget(self.slicePlotCheckboxWidget)
        self.slicePlotAxisLayout.addWidget(self.slicePlotSliceWidthWidget)
        # self.slicePlotXAxisCombo.currentIndexChanged.connect(self.plotDataSlice)
        self.slicePlotSliceWidthWidget.valueChanged.connect(self.changeSliceLength)
        # self.beamPlotXAxisCombo.setCurrentIndex(2)
        # self.beamPlotYAxisCombo.setCurrentIndex(5)
        self.slicePlotLayout.addWidget(self.slicePlotAxisWidget)
        self.slicePlotLayout.addWidget(self.slicePlotWidgetGraphicsLayout)

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
        self.tabWidget.addTab(self.slicePlotWidget,'Slice Beam Plots')
        self.tabWidget.currentChanged.connect(self.changeTab)
        self.layout.addWidget(self.folderBeamWidget)
        self.layout.addWidget(self.tabWidget)

        self.plotType = 'Twiss'
        self.changeDirectory(self.directory)

    def findFirstEmptyColumnInGraphicsLayout(self):
            rowsfilled =  self.slicePlotWidgetGraphicsLayout.ci.rows.get(0, {}).keys()
            for i in range(49):
                if not i in rowsfilled:
                    return i

    def changeTab(self, i):
        if self.tabWidget.tabText(i) == 'Beam Plots':
            self.plotType = 'Beam'
            self.beamWidget.setVisible(True)
        elif self.tabWidget.tabText(i) == 'Slice Beam Plots':
            self.plotType = 'Slice'
            self.beamWidget.setVisible(True)
        else:
            self.plotType = 'Twiss'
            self.beamWidget.setVisible(False)
        self.loadDataFile()

    def changeDirectory(self, directory=None):
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
        self.loadDataFile()

    def getScreenFiles(self):
        self.screenpositions = {}
        files = glob.glob(self.directory+'/*.????.???')
        filenames = ['.'.join(os.path.basename(f).split('.')[:2]) for f in files]
        runnumber = [os.path.basename(f).split('.')[-1] for f in files]
        for f, r in list(set(zip(filenames, runnumber))):
            files = glob.glob(self.directory+'/'+f+'.????.???')
            screenpositions = [re.search('\d\d\d\d', s).group(0) for s in files]
            self.screenpositions[f] = {'screenpositions': sorted(screenpositions), 'run': r}

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
        if self.plotType == 'Twiss':
            files = sorted(glob.glob(self.directory+"/*Xemit*"))
            self.twiss.read_astra_emit_files(files)
            self.plotDataTwiss()
        elif self.plotType == 'Beam' or self.plotType == 'Slice':
            if hasattr(self,'beamFileName') and os.path.isfile(self.directory+'/'+self.beamFileName):
                self.beam.read_astra_beam_file(self.directory+'/'+self.beamFileName)
                if self.plotType == 'Beam':
                    self.plotDataBeam()
                else:
                    self.beam.bin_time()
                    self.plotDataSlice()

    def plotDataTwiss(self):
        for entry in self.twissplotLayout:
            if entry == 'next_row':
                pass
            else:
                x = self.twiss['z']
                y = self.twiss[entry['name']]*entry['scale']
                self.twissPlots[entry['name']].setData(x=x, y=y, pen=mkPen('b', width=3))

    def plotDataBeam(self):
        self.beamPlot.setData(x=getattr(self.beam, str(self.beamPlotXAxisCombo.currentText())), y=getattr(self.beam, str(self.beamPlotYAxisCombo.currentText())), pen=None, symbol='+')

    def changeSliceLength(self):
        self.beam.slice_length = self.slicePlotSliceWidthWidget.value()*1e-15
        self.beam.bin_time()
        self.plotDataSlice()

    def plotDataSlice(self):
        for param in self.sliceParams:
            if self.slicePlotCheckbox[param['name']].isChecked():
                x = self.beam.slice_bins
                self.slicePlot.setRange(xRange=[min(x),max(x)])
                # self.plot.setRange(xRange=[-0.5,1.5])
                y = getattr(self.beam, param['name'])
                self.curve[param['name']].setData(x=x, y=y)
                self.sliceaxis[param['name']][0].setVisible(True)
                # currentrange = self.sliceaxis[param['name']][0].range
                # print 'currentrange = ', currentrange
                # self.sliceaxis[param['name']][0].setRange(0, currentrange[1])
            else:
                # pass
                self.curve[param['name']].setData(x=[], y=[])
                self.sliceaxis[param['name']][0].setVisible(False)
            self.sliceaxis[param['name']][1].autoRange()
            currentrange = self.sliceaxis[param['name']][1].viewRange()
            self.sliceaxis[param['name']][1].setYRange(0, currentrange[1][1])

    def scaleLattice(self, vb, range):
        yrange = range[1]
        scaleY = 0.05*abs(yrange[1] - yrange[0])
        rect = QRectF(0, yrange[0] + 2*scaleY, self.twiss['z'][-1], 4*scaleY)
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
