import sys, os, time, math, datetime, copy, re,  h5py
from copy import copy
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
from SimulationFramework.Modules.multiAxisPlot import multiAxisPlotWidget
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

class slicePlotWidget(multiAxisPlotWidget):
    ''' QWidget containing pyqtgraph plot showing slice beam parameters '''

    plotParams = [
        {'label': 'Horizontal Emittance (normalised)', 'quantity': 'slice_normalized_horizontal_emittance', 'units': 'm-rad', 'name': '&epsilon;<sub>n,x</sub>'},
        {'label': 'Vertical Emittance (normalised)', 'quantity': 'slice_normalized_vertical_emittance', 'units': 'm-rad', 'name': '&epsilon;<sub>n,y</sub>'},
        {'label': 'Current', 'quantity': 'slice_peak_current', 'units': 'A', 'text': 'I', 'name': 'I'},
        {'label': 'Relative Momentum Spread', 'quantity': 'slice_relative_momentum_spread', 'units': '%', 'name': '&sigma;<sub>cp</sub>/p'},
        {'label': 'Horizontal Beta Function', 'quantity': 'slice_beta_x', 'units': 'm', 'name': '&beta;<sub>x</sub>'},
        {'label': 'Vertical Beta Function', 'quantity': 'slice_beta_y', 'units': 'm', 'name': '&beta;<sub>y</sub>'},
    ]

    def __init__(self, **kwargs):
        super(slicePlotWidget, self).__init__(**kwargs)
        self.beams = {}
        self.slicePlotSliceWidthWidget = QSpinBox()
        self.slicePlotSliceWidthWidget.setMaximum(500)
        self.slicePlotSliceWidthWidget.setValue(100)
        self.slicePlotSliceWidthWidget.setSingleStep(10)
        self.slicePlotSliceWidthWidget.setSuffix(" slices")
        self.slicePlotSliceWidthWidget.setSpecialValueText('Automatic')
        self.multiaxisPlotAxisLayout.addWidget(self.slicePlotSliceWidthWidget)
        self.slicePlotSliceWidthWidget.valueChanged.connect(self.changeSliceLengths)

    def addsliceDataFiles(self, dicts):
        '''
            add multiple data dictionaries

            Keyword arguments:
            dicts -- dictionary containing directory definitions:
                [
                    {'directory': <dir location>,           'filename': [<list of HDF5 beam files>]},
                    ...
                ]
        '''
        for d in dicts:
            self.addsliceDataFile(**d)

    def addsliceDataObject(self, beamobject, name):
        '''
            addsliceDirectory - read the data files in a directory and add a plotItem to the relevant slicePlotItems

            Keyword arguments:
            beamobject -- beam object
            name -- key index name
        '''
        ''' load the data files into the slice dictionary '''
        if str(type(beamobject)) == "<class 'SimulationFramework.Modules.read_beam_file.beam'>":
            self.beams[name] = beamobject
            beamobject.bin_time()
            self.curves[name] = {}
            for n, param in enumerate(self.plotParams):
                label = param['label']
                color = self.colors[n % len(self.colors)]
                pen = pg.mkPen(color=color, style=self.styles[self.plotColor % len(self.styles)], width=3)
                exponent = np.floor(np.log10(np.abs(beamobject.slice_length)))
                x = 10**(12) * np.array((beamobject.slice_bins - np.mean(beamobject.slice_bins)))
                # self.multiaxisPlot.setRange(xRange=[min(x),max(x)])
                y = getattr(beamobject, param['quantity'])
                self.addCurve(x, y, name, label, pen)
            self.plotColor += 1
        self.updateMultiAxisPlot()

    def addsliceDataFile(self, directory, filename=None):
        '''
            addsliceDirectory - read the data files in a directory and call addsliceDataObject on the beam object

            Keyword arguments:
            directory -- location of beam file(s)
            filenames -- HDF5 filenames of beam file(s)
        '''
        if not isinstance(filename, (list, tuple)):
            filename = [filename]
        for f in filename:
            datafile = directory + '/' + f if f is not None else directory
            if os.path.isfile(datafile):
                beam = raf.beam()
                beam.read_HDF5_beam_file(datafile)
                self.addsliceDataObject(beam, datafile)

    def changeSliceLengths(self):
        ''' change the time slice length for all beam objects '''
        for d in self.beams:
            self.changeSliceLength(d)
        self.updateMultiAxisPlot()

    def changeSliceLength(self, datafile):
        ''' change the time slice length for a beam data object and update the plot '''
        beam = self.beams[datafile]
        beam.slices = self.slicePlotSliceWidthWidget.value()
        beam.bin_time()
        for n, param in enumerate(self.plotParams):
            label = param['label']
            exponent = np.floor(np.log10(np.abs(beam.slice_length)))
            x = 10**(12) * np.array((beam.slice_bins - np.mean(beam.slice_bins)))
            self.multiaxisPlot.setRange(xRange=[min(x),max(x)])
            y = getattr(beam, param['quantity'])
            self.curves[datafile][label].setData(x=x, y=y)

def main():
    app = QApplication(sys.argv)
    pg.setConfigOptions(antialias=True)
    pg.setConfigOption('background', 'w')
    pg.setConfigOption('foreground', 'k')
    ex = mainWindow()
    ex.show()
    ex.slicePlot.addsliceDataFiles([
    {'directory': 'OnlineModel_test_data/basefiles_4_250pC', 'filename': 'CLA-S02-APER-01.hdf5'},
    {'directory': 'OnlineModel_test_data/test_4', 'filename': ['CLA-L02-APER.hdf5','CLA-S04-APER-01.hdf5']}])
    ex.slicePlot.addsliceDataFile('OnlineModel_test_data/test_4/CLA-S03-APER.hdf5')
    ex.slicePlot.highlightPlot('OnlineModel_test_data/test_4/CLA-S03-APER.hdf5')
    sys.exit(app.exec_())

if __name__ == '__main__':
   main()
