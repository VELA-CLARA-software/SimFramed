import sys, os, time, math, datetime, copy, re,  h5py
from collections import OrderedDict
import glob
from PyQt5.QtCore import * #QObject, pyqtSignal, QThread, QTimer, QRectF, Qt
from PyQt5.QtGui import * #QApplication, QMainWindow, QWidget, QPushButton, QVBoxLayout, QHBoxLayout, QComboBox, QTabWidget, QLineEdit, QFileDialog, QLabel, QAction, QPixmap, qApp, QStyle, QGroupBox, QSpinBox
from PyQt5.QtWidgets import *
from pyqtgraph import LegendItem, mkPen, mkBrush, LabelItem, TableWidget, GraphicsLayoutWidget, setConfigOption, \
setConfigOptions, InfiniteLine, ImageItem, GraphicsView, GraphicsLayout, AxisItem, ViewBox, PlotDataItem, colorStr, mkColor, ImageView, PlotItem
from pyqtgraph.graphicsItems.LegendItem import ItemSample
import argparse
import imageio
import numpy as np
sys.path.append(os.path.abspath(os.path.realpath(__file__)+'/../../../'))
print (sys.path)
import SimulationFramework.Modules.read_beam_file as raf
import SimulationFramework.Modules.read_twiss_file as rtf
sys.path.append(os.path.realpath(__file__)+'/../../../../')
# import Software.Widgets.loggerWidget.loggerWidget as lw
# import logging
# logger = logging.getLogger(__name__)

from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib.colors import LogNorm

# try:
#      = QString.fromUtf8
# except AttributeError:
#     def (s):
#         return s

try:
    _encoding = QApplication.UnicodeUTF8
    def _translate(context, text, disambig):
        return QApplication.translate(context, text, disambig, _encoding)
except AttributeError:
    def _translate(context, text, disambig):
        return QApplication.translate(context, text, disambig)


class MyMplCanvas(FigureCanvas):
    """Ultimately, this is a QWidget (as well as a FigureCanvasAgg, etc.)."""

    def __init__(self, parent=None, width=5, height=5, dpi=50):
        fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes = fig.add_subplot(111)

        self.compute_initial_figure()

        FigureCanvas.__init__(self, fig)
        self.setParent(parent)
        # FigureCanvas.setSizePolicy(self,
        #                            QSizePolicy.Expanding,
        #                            QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)

    def compute_initial_figure(self):
        pass

class MyStaticMplCanvas(MyMplCanvas):
    """Simple canvas with a sine plot."""

    def plothist(self, x, y, nbins, xLabel=None, yLabel=None):
        self.axes.cla()
        # self.axes.set_aspect('equal', 'datalim')
        self.axes.hist2d(x, y, bins=nbins, norm=LogNorm())
        if xLabel is not None:
            self.setXlabel(xLabel)
        if yLabel is not None:
            self.setYlabel(yLabel)
        self.draw()

    def setXlabel(self, lbl):
        self.axes.set_xlabel(lbl)

    def setYlabel(self, lbl):
        self.axes.set_ylabel(lbl)

parser = argparse.ArgumentParser(description='Plot ASTRA Data Files')
parser.add_argument('-d', '--directory', default='.')

def rainbow():
    array = np.empty((256, 3))
    abytes = np.arange(0, 1, 0.00390625)
    array[:, 0] = np.abs(2 * abytes - 0.5) * 255
    array[:, 1] = np.sin(abytes * np.pi) * 255
    array[:, 2] = np.cos(abytes * np.pi / 2) * 255
    return array

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

        self.latticePlotData = imageio.imread(os.path.dirname( os.path.abspath(__file__))+'/lattice_plot.png')
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
                # latticePlot.setOpacity(0.5)
                self.twissPlots[entry['name']] = p.plot(pen=mkPen('b', width=3))
                self.latticePlots[p.vb] = latticePlot
                p.vb.sigRangeChanged.connect(self.scaleLattice)

        ''' beamPlotWidget '''
        self.beamPlotWidget = QWidget()
        self.beamPlotLayout = QVBoxLayout()
        self.beamPlotWidget.setLayout(self.beamPlotLayout)

        #
        # self.beamPlotChoice =
        #

        self.beamPlotAxisWidget = QWidget()
        self.beamPlotAxisWidget.setMaximumHeight(100)
        Form = self.beamPlotAxisWidget

        self.beamPlotXAxisDict = OrderedDict()
        self.beamPlotXAxisDict['x'] = {'scale':1e3, 'axis': 'x [mm]'}
        self.beamPlotXAxisDict['y'] = {'scale':1e3, 'axis': 'y [mm]'}
        self.beamPlotXAxisDict['z'] = {'scale':1e6, 'axis': 'z [micron]', 'norm': True}
        self.beamPlotXAxisDict['t'] = {'scale':1e12, 'axis': 't [ps]', 'norm': True}
        self.beamPlotXAxisDict['cpx'] = {'scale':1e3, 'axis': 'cpx [keV]'}
        self.beamPlotXAxisDict['cpy'] = {'scale':1e3, 'axis': 'cpy [keV]'}
        self.beamPlotXAxisDict['BetaGamma']= {'scale':0.511 , 'axis': 'cp [MeV]'}

        # Form.setObjectName(("Form"))
        # Form.resize(874, 212)
        sizePolicy = QSizePolicy(QSizePolicy.Preferred, QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(Form.sizePolicy().hasHeightForWidth())
        Form.setSizePolicy(sizePolicy)
        self.horizontalLayout = QHBoxLayout(Form)
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.beamPlotXAxisCombo = QComboBox(Form)
        self.beamPlotXAxisCombo.addItems(list(self.beamPlotXAxisDict.keys()))
        self.beamPlotXAxisCombo.setCurrentIndex(2)
        self.horizontalLayout.addWidget(self.beamPlotXAxisCombo)
        self.beamPlotYAxisCombo = QComboBox(Form)
        self.beamPlotYAxisCombo.addItems(list(self.beamPlotXAxisDict.keys()))
        self.beamPlotYAxisCombo.setCurrentIndex(6)
        self.horizontalLayout.addWidget(self.beamPlotYAxisCombo)
        self.groupBox = QGroupBox(Form)
        self.groupBox.setObjectName("groupBox")
        self.formLayout = QFormLayout(self.groupBox)
        self.formLayout.setFieldGrowthPolicy(QFormLayout.AllNonFixedFieldsGrow)
        self.formLayout.setObjectName("formLayout")
        self.chooseSliceWidth = QRadioButton(self.groupBox)
        self.chooseSliceWidth.setObjectName(("chooseSliceWidth"))
        self.formLayout.setWidget(0, QFormLayout.LabelRole, self.chooseSliceWidth)
        self.sliceWidth = QDoubleSpinBox(self.groupBox)
        self.sliceWidth.setDecimals(6)
        self.sliceWidth.setObjectName("sliceWidth")
        self.formLayout.setWidget(0, QFormLayout.FieldRole, self.sliceWidth)
        self.chooseSliceNumber = QRadioButton(self.groupBox)
        self.chooseSliceNumber.setChecked(True)
        self.chooseSliceNumber.setObjectName("chooseSliceNumber")
        self.formLayout.setWidget(1, QFormLayout.LabelRole, self.chooseSliceNumber)
        self.sliceNumber = QSpinBox(self.groupBox)
        self.sliceNumber.setObjectName(("sliceNumber"))
        self.formLayout.setWidget(1, QFormLayout.FieldRole, self.sliceNumber)
        self.horizontalLayout.addWidget(self.groupBox)
        self.chooseSliceWidth.setText(_translate("Form", "Slice Width", None))
        self.chooseSliceNumber.setText(_translate("Form", "Number of Slices", None))
        self.sliceWidth.setRange(1e-6,1)
        self.sliceWidth.setSingleStep(0.00001)
        self.histogramWidth = 0.0005
        self.sliceWidth.setValue(self.histogramWidth)
        # self.sliceNumber = QSpinBox()
        self.sliceNumber.setRange(10,10000)
        self.sliceNumber.setSingleStep(10)
        self.histogramBins = 100
        self.sliceNumber.setValue(self.histogramBins)

        # self.beamPlotXAxisCombo = QComboBox()

        # self.beamPlotAxisLayout = QHBoxLayout()
        # self.beamPlotAxisWidget.setLayout(self.beamPlotAxisLayout)

        # self.beamPlotAxisLayout.addWidget(self.beamPlotXAxisCombo)
        # self.beamPlotAxisLayout.addWidget(self.beamPlotYAxisCombo)
        # self.beamPlotAxisLayout.addWidget(self.beamPlotNumberBins)
        self.beamPlotXAxisCombo.currentIndexChanged.connect(self.plotDataBeam)
        self.beamPlotYAxisCombo.currentIndexChanged.connect(self.plotDataBeam)
        self.chooseSliceNumber.toggled.connect(self.plotDataBeam)
        self.sliceNumber.valueChanged.connect(self.plotDataBeam)
        self.sliceWidth.valueChanged.connect(self.plotDataBeam)
            # self.beamPlotXAxisCombo.setCurrentIndex(2)
            # self.beamPlotYAxisCombo.setCurrentIndex(5)
        self.canvasWidget = QWidget()
        l = QVBoxLayout(self.canvasWidget)
        self.sc = MyStaticMplCanvas(self.canvasWidget, width=1, height=1, dpi=150)
        l.addWidget(self.sc)
        self.beamPlotLayout.addWidget(self.beamPlotAxisWidget)
        self.beamPlotLayout.addWidget(self.canvasWidget)

        ''' slicePlotWidget '''
        self.sliceParams = [
        {'name': 'slice_normalized_horizontal_emittance', 'units': 'm-rad', 'text': 'enx'},
        {'name': 'slice_normalized_vertical_emittance', 'units': 'm-rad', 'text': 'eny'},
        {'name': 'slice_peak_current', 'units': 'A', 'text': 'PeakI'},
        {'name': 'slice_relative_momentum_spread', 'units': '%', 'text': 'sigma-p'},
        {'name': 'slice_beta_x', 'units': 'm', 'text': 'beta_x'},
        {'name': 'slice_beta_y', 'units': 'm', 'text': 'beta_y'},
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
        colors = ['b','r','g','k','y','m','c']
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
        self.slicePlotSliceWidthWidget.setMaximum(500)
        self.slicePlotSliceWidthWidget.setValue(100)
        self.slicePlotSliceWidthWidget.setSingleStep(10)
        self.slicePlotSliceWidthWidget.setSuffix(" slices")
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
        self.reloadButton.clicked.connect(lambda: self.changeDirectory(self.directory))
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
        # self.log = lw.loggerWidget()
        # self.log.addLogger(logger)
        # sys.stdout = lw.redirectLogger(self.log, 'stdout')
        # self.tabWidget.addTab(self.log,'Log')
        self.tabWidget.currentChanged.connect(self.changeTab)
        self.layout.addWidget(self.folderBeamWidget)
        self.layout.addWidget(self.tabWidget)

        self.plotType = 'Twiss'
        self.changeDirectory(self.directory)

    def findFirstEmptyColumnInGraphicsLayout(self):
            rowsfilled =  list(self.slicePlotWidgetGraphicsLayout.ci.rows.get(0, {}).keys())
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
        self.screenSelector.currentIndexChanged.disconnect(self.changeScreen)
        self.fileSelector.currentIndexChanged.disconnect(self.updateScreenCombo)
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
        self.screenSelector.currentIndexChanged.connect(self.changeScreen)
        self.fileSelector.currentIndexChanged.connect(self.updateScreenCombo)
        self.changeScreen(self.screenSelector.currentIndex())

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
        # print 'filenames = ', filenames
        runnumber = ['001' for f in filenames]
        # print 'runnumber = ', runnumber
        for f, r in list(set(zip(filenames, runnumber))):
            files = glob.glob(self.directory+'/'+f+'*.hdf5')
            screennames = sorted([os.path.basename(s).split('.')[0] for s in files], key=lambda x: self.getSposition(x))
            screenpositions = [self.getSposition(s) for s in screennames]
            # print 'screenpositions = ', screenpositions
            self.screenpositions[f] = {'screenpositions': screenpositions, 'screennames': screennames, 'run': r}

    def updateFileCombo(self):
        self.fileSelector.clear()
        i = -1
        screenfirstpos = []
        for f in self.screenpositions:
            if len(self.screenpositions[f]['screenpositions']) > 0:
                screenfirstpos.append([f, min(self.screenpositions[f]['screenpositions'])])
        screenfirstpos = np.array(screenfirstpos)
        # print 'screenfirstpos = ', screenfirstpos
        sortedscreennames = sorted(screenfirstpos, key=lambda x: float(x[1]))
        print('sortedscreennames = ', sortedscreennames)
        for f in sortedscreennames:
            self.fileSelector.addItem(f[0])
            i += 1
            if f[0] == self.currentFileText:
                self.fileSelector.setCurrentIndex(i)

    def changeScreen(self, i):
        run = self.screenpositions[str(self.fileSelector.currentText())]['run']
        self.beamFileName = str(self.screenpositions[str(self.fileSelector.currentText())]['screennames'][self.screenSelector.currentIndex()])+'.hdf5'
        print('beamFileName = ', self.beamFileName)
        self.loadDataFile()

    def updateScreenCombo(self):
        self.screenSelector.clear()
        i = -1
        for i, s in enumerate(self.screenpositions[str(self.fileSelector.currentText())]['screennames']):
            self.screenSelector.addItem(str(s)+' ('+str(self.screenpositions[str(self.fileSelector.currentText())]['screenpositions'][i])+')')
            i += 1
            if s == self.currentScreenText:
                self.screenSelector.setCurrentIndex(i)

    def loadDataFile(self):
        if self.plotType == 'Twiss':
            files = sorted(glob.glob(self.directory+"/*Xemit*"))
            self.twiss.read_astra_emit_files(files)
            files = sorted(glob.glob(self.directory+"/*.flr"))
            self.twiss.read_elegant_twiss_files(files, reset=False)
            self.plotDataTwiss()
        elif self.plotType == 'Beam' or self.plotType == 'Slice':
            if hasattr(self,'beamFileName') and os.path.isfile(self.directory+'/'+self.beamFileName):
                # starttime = time.time()
                self.beam.read_HDF5_beam_file(self.directory+'/'+self.beamFileName)
                # print 'reading file took ', time.time()-starttime, 's'
                # print 'Read file: ', self.beamFileName
                if self.plotType == 'Beam':
                    self.plotDataBeam()
                else:
                    self.changeSliceLength()
                    # self.plotDataSlice()

    def plotDataTwiss(self):
        for entry in self.twissplotLayout:
            if entry == 'next_row':
                pass
            else:
                x = self.twiss['z']
                y = self.twiss[entry['name']]*entry['scale']
                xy = np.transpose(np.array([x,y]))
                x, y = np.transpose(xy[np.argsort(xy[:,0])])
                self.twissPlots[entry['name']].setData(x=x, y=y, pen=mkPen('b', width=3))

    def plotDataBeam(self):
        xdict = self.beamPlotXAxisDict[str(self.beamPlotXAxisCombo.currentText())]
        ydict = self.beamPlotXAxisDict[str(self.beamPlotYAxisCombo.currentText())]
        x = xdict['scale'] * getattr(self.beam, str(self.beamPlotXAxisCombo.currentText()))
        if 'norm' in xdict and xdict['norm'] == True:
            x = x - np.mean(x)
        y = ydict['scale'] * getattr(self.beam, str(self.beamPlotYAxisCombo.currentText()))
        if 'norm' in ydict and ydict['norm'] == True:
            y = y - np.mean(y)
        if self.chooseSliceWidth.isChecked():
            slices = np.arange(min(x), max(x) + self.sliceWidth.value(), self.sliceWidth.value())
            self.sliceNumber.setValue(len(slices))
            self.sc.plothist(x,y, [slices, len(slices)], xLabel=xdict['axis'], yLabel=ydict['axis'])
        else:
            self.sliceWidth.setValue(int(1e6*((max(x) - min(x)) / self.sliceNumber.value()))/1e6)
            self.sc.plothist(x,y, self.sliceNumber.value(), xLabel=xdict['axis'], yLabel=ydict['axis'])

    def changeSliceLength(self):
        self.beam.slices = self.slicePlotSliceWidthWidget.value()
        self.beam.bin_time()
        self.plotDataSlice()

    def plotDataSlice(self):
        for param in self.sliceParams:
            if self.slicePlotCheckbox[param['name']].isChecked():
                exponent = np.floor(np.log10(np.abs(self.beam.slice_length)))
                x = 10**(12) * np.array((self.beam.slice_bins - np.mean(self.beam.slice_bins)))
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
        rect = QRectF(0, yrange[0] + 2*scaleY, 49.2778, 4*scaleY)
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
