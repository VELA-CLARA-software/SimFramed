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
import SimulationFramework.Modules.read_beam_file as raf
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
        self.slicePlot = slicePlotWidget()

        self.layout.addWidget(self.slicePlot)

        self.setCentralWidget(self.centralWidget)

        self.setWindowTitle("ASTRA Data Plotter")
        menubar = self.menuBar()
        fileMenu = menubar.addMenu('&File')

        exitAction = QAction('&Exit', self)
        exitAction.setShortcut('Ctrl+Q')
        exitAction.setStatusTip('Exit application')
        exitAction.triggered.connect(self.close)
        fileMenu.addAction(exitAction)

class slicePlotWidget(QWidget):
    # Styles for the plot lines
    colors = [QColor('#F5973A'),QColor('#A95AA1'),QColor('#85C059'),QColor('#0F2080'),QColor('#BDB8AD'), 'r']
    styles = [Qt.SolidLine, Qt.DashLine, Qt.DotLine, Qt.DashDotLine, Qt.DashDotDotLine]

    sliceParams = [
        {'label': 'Horizontal Emittance (normalised)', 'quantity': 'slice_normalized_horizontal_emittance', 'units': 'm-rad', 'name': '&epsilon;<sub>n,x</sub>'},
        {'label': 'Vertical Emittance (normalised)', 'quantity': 'slice_normalized_vertical_emittance', 'units': 'm-rad', 'name': '&epsilon;<sub>n,y</sub>'},
        {'label': 'Current', 'quantity': 'slice_peak_current', 'units': 'A', 'text': 'I', 'name': 'I'},
        {'label': 'Relative Momentum Spread', 'quantity': 'slice_relative_momentum_spread', 'units': '%', 'name': '&sigma;<sub>cp</sub>/p'},
        {'label': 'Horizontal Beta Function', 'quantity': 'slice_beta_x', 'units': 'm', 'name': '&beta;<sub>x</sub>'},
        {'label': 'Vertical Beta Function', 'quantity': 'slice_beta_y', 'units': 'm', 'name': '&beta;<sub>y</sub>'},
    ]

    def __init__(self, **kwargs):
        super(slicePlotWidget, self).__init__(**kwargs)
        ''' These are for reading data files from ASTRA and Elegant '''
        self.beams = {}
        self.twiss = rtf.twiss()

        self.slicePlotWidget = QWidget()
        self.slicePlotLayout = QVBoxLayout()
        self.slicePlotWidget.setLayout(self.slicePlotLayout)
        self.slicePlotWidgetGraphicsLayout = pg.GraphicsLayoutWidget()
        self.slicePlotCheckbox = {}
        self.viewboxes = {}
        self.curves = {}
        self.sliceaxis = {}
        self.slicePlotCheckboxWidget = QWidget()
        self.slicePlotCheckboxLayout = QVBoxLayout()
        self.slicePlotCheckboxWidget.setLayout(self.slicePlotCheckboxLayout)
        self.slicePlot = self.slicePlotWidgetGraphicsLayout.addPlot(title='Slice',row=0,col=50)
        self.slicePlot.showAxis('left', False)
        self.slicePlot.showGrid(x=True, y=True)
        i = 0
        for param in self.sliceParams:
            axis = pg.AxisItem("left")
            labelStyle = {'color': '#'+pg.colorStr(pg.mkColor(self.colors[i]))[0:-2]}
            axis.setLabel(text=param['name'], units=param['units'], **labelStyle)
            i += 1
            viewbox = pg.ViewBox()
            self.viewboxes[param['label']] = viewbox
            axis.linkToView(viewbox)
            viewbox.setXLink(self.slicePlot.vb)
            self.sliceaxis[param['label']] = [axis, viewbox]
            self.curves[param['label']] = {}
            col = i
            self.slicePlotWidgetGraphicsLayout.ci.addItem(axis, row = 0, col = col,  rowspan=1, colspan=1)
            self.slicePlotWidgetGraphicsLayout.ci.addItem(viewbox, row=0, col=50)
            viewbox.setLimits(yMin=0)
            self.slicePlotCheckbox[param['label']] = QCheckBox(param['label'])
            self.slicePlotCheckbox[param['label']].setChecked(True)
            self.slicePlotCheckboxLayout.addWidget(self.slicePlotCheckbox[param['label']])
            self.slicePlotCheckbox[param['label']].stateChanged.connect(self.updateSlicePlot)
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
        self.slicePlotSliceWidthWidget.valueChanged.connect(self.changeSliceLengths)
        self.slicePlotLayout.addWidget(self.slicePlotAxisWidget)
        self.slicePlotLayout.addWidget(self.slicePlotWidgetGraphicsLayout)

        self.layout = QVBoxLayout()
        self.setLayout(self.layout)

        self.layout.addWidget(self.slicePlotWidget)

        ''' used for style cycling '''
        self.plotColor = 0

    def addsliceDataFiles(self, dicts):
        for d in dicts:
            self.addsliceDataFile(**d)

    def addsliceDataObject(self, beamobject, datafile):
        ''' addsliceDirectory - read the data files in a directory and add a plotItem to the relevant slicePlotItems '''
        ''' load the data files into the slice dictionary '''
        if str(type(beamobject)) == "<class 'SimulationFramework.Modules.read_beam_file.beam'>":
            self.beams[datafile] = beamobject
            beamobject.bin_time()
            for n, param in enumerate(self.sliceParams):
                color = self.colors[n]
                pen = pg.mkPen(color=color, style=self.styles[self.plotColor % len(self.styles)])
                self.curves[param['label']][datafile] = pg.PlotDataItem()
                self.viewboxes[param['label']].addItem(self.curves[param['label']][datafile])
                exponent = np.floor(np.log10(np.abs(beamobject.slice_length)))
                x = 10**(12) * np.array((beamobject.slice_bins - np.mean(beamobject.slice_bins)))
                self.slicePlot.setRange(xRange=[min(x),max(x)])
                y = getattr(beamobject, param['quantity'])
                self.curves[param['label']][datafile].setData(x=x, y=y, pen=pen)
            self.plotColor += 1
        self.updateSlicePlot()

    def addsliceDataFile(self, directory, filename):
        ''' addsliceDirectory - read the data files in a directory and add a plotItem to the relevant slicePlotItems '''
        ''' load the data files into the slice dictionary '''
        datafile = directory + '/' + filename
        if os.path.isfile(datafile):
            beam = raf.beam()
            beam.read_HDF5_beam_file(datafile)
            self.addsliceDataObject(beam, datafile)

    def changeSliceLengths(self):
        for d in self.beams:
            self.changeSliceLength(d)

    def changeSliceLength(self, datafile):
        beam = self.beams[datafile]
        beam.slices = self.slicePlotSliceWidthWidget.value()
        beam.bin_time()
        for n, param in enumerate(self.sliceParams):
            exponent = np.floor(np.log10(np.abs(beam.slice_length)))
            x = 10**(12) * np.array((beam.slice_bins - np.mean(beam.slice_bins)))
            self.slicePlot.setRange(xRange=[min(x),max(x)])
            y = getattr(beam, param['quantity'])
            self.curves[param['label']][datafile].setData(x=x, y=y)

    def updateSlicePlot(self):
        for n, param in enumerate(self.sliceParams):
            if self.slicePlotCheckbox[param['label']].isChecked():
                for d in self.curves[param['label']]:
                    self.curves[param['label']][d].setVisible(True)
                self.sliceaxis[param['label']][0].setVisible(True)
            else:
                for d in self.curves[param['label']]:
                    self.curves[param['label']][d].setVisible(False)
                self.sliceaxis[param['label']][0].setVisible(False)
            self.sliceaxis[param['label']][1].autoRange()
            currentrange = self.sliceaxis[param['label']][1].viewRange()
            self.sliceaxis[param['label']][1].setYRange(0, currentrange[1][1])

    def removePlot(self, directory):
        ''' finds all slice plots based on a directory name, and removes them '''
        if not isinstance(directory, (list, tuple)):
            directory = [directory]
        indexname = directory[-1]
        if directory in self.slicePlotItems:
            for entry in self.sliceplotLayout:
                if entry == 'next_row':
                    pass
                else:
                    self.slicePlotWidgets[entry['label']].removeItem(self.slicePlotItems[indexname][entry['label']])

pg.setConfigOptions(antialias=True)
pg.setConfigOption('background', 'w')
pg.setConfigOption('foreground', 'k')
def main():
    app = QApplication(sys.argv)
    pg.setConfigOptions(antialias=True)
    pg.setConfigOption('background', 'w')
    pg.setConfigOption('foreground', 'k')
    ex = mainWindow()
    ex.show()
    ex.slicePlot.addsliceDataFiles([
    {'directory': 'OnlineModel_test_data/basefiles_4_250pC', 'filename': 'CLA-S02-APER-01.hdf5'},
    {'directory': 'OnlineModel_test_data/test_4', 'filename': 'CLA-L02-APER.hdf5'}])
    sys.exit(app.exec_())

if __name__ == '__main__':
   main()
