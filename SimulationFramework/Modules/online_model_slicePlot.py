import sys, os, time, math, datetime, copy, re,  h5py
from collections import OrderedDict
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
from SimulationFramework.Modules.multiPlot import multiPlotWidget

class slicePlotter(QMainWindow):
    def __init__(self):
        super(slicePlotter, self).__init__()
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

class slicePlotWidget(multiPlotWidget):
    # Layout oder for the individual Tiwss plot items
    plotParams = [
        {'label': 'Slice Horizontal Emittance (normalised)', 'quantity': 'slice_normalized_horizontal_emittance', 'units': 'm-rad', 'name': '&epsilon;<sub>n,x</sub>', 'range': [0,5e-6]},
        {'label': 'Slice Vertical Emittance (normalised)', 'quantity': 'slice_normalized_vertical_emittance', 'units': 'm-rad', 'name': '&epsilon;<sub>n,y</sub>', 'range': [0,5e-6]},
        'next_row',
        {'label': 'Current', 'quantity': 'slice_peak_current', 'units': 'A', 'text': 'I', 'name': 'I', 'range': [0,100]},
        {'label': 'Slice Relative Momentum Spread', 'quantity': 'slice_relative_momentum_spread', 'units': '%', 'name': '&sigma;<sub>cp</sub>/p', 'range': [0,0.2]},
        'next_row',
        {'label': 'Slice Horizontal Beta Function', 'quantity': 'slice_beta_x', 'units': 'm', 'name': '&beta;<sub>x</sub>', 'range': [0,250]},
        {'label': 'Slice Vertical Beta Function', 'quantity': 'slice_beta_y', 'units': 'm', 'name': '&beta;<sub>y</sub>', 'range': [0,250]},
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
        self.slicePlotSliceWidthWidget.setMaximumWidth(150)
        # self.multiaxisPlotAxisLayout.addWidget(self.slicePlotSliceWidthWidget)
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

    def addsliceDataObject(self, beamobject, id, color=None):
        '''
            addsliceDirectory - read the data files in a directory and add a plotItem to the relevant slicePlotItems

            Keyword arguments:
            beamobject -- beam object
            name -- key index name
        '''
        ''' load the data files into the slice dictionary '''
        if str(type(beamobject)) == "<class 'SimulationFramework.Modules.read_beam_file.beam'>":
            self.beams[id] = beamobject
            self.beams[id].slices = self.slicePlotSliceWidthWidget.value()
            beamobject.bin_time()
            if color is None:
                color = self.colors[self.plotColor % len(self.colors)]
                if id not in self.curves:
                    self.plotColor += 1
            pen = pg.mkPen(color, width=3, style=self.styles[int(self.plotColor % len(self.styles))])
            if id not in self.curves:
                self.curves[id] = {}
            for n, param in enumerate(self.plotParams):
                if not param == 'next_row':
                    label = param['label']
                    # color = self.colors[n]
                    # pen = pg.mkPen(color=color, style=self.styles[self.plotColor % len(self.styles)], width=3)
                    exponent = np.floor(np.log10(np.abs(beamobject.slice_length)))
                    x = 10**(12) * np.array((beamobject.slice_bins - np.mean(beamobject.slice_bins)))
                    # self.multiPlot.setRange(xRange=[min(x),max(x)])
                    y = getattr(beamobject, param['quantity'])
                    # print('#################################################################')
                    # print(name, name in self.curves, label, self.curves)#, label, label in self.curves[name], self.curves[name])
                    # print('#################################################################')
                    if id in self.curves and label in self.curves[id]:
                        # print('Updating curve: ', name, label)
                        self.updateCurve(x, y, id, label)
                    else:
                        # print('ADDING curve: ', name, label)
                        self.addCurve(x, y, id, label, pen)


    def addsliceDataFile(self, directory, filename=None, color=None, id=None):
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
            if id is None:
                id = datafile
            if not datafile in self.curves and os.path.isfile(datafile):
                beam = raf.beam()
                # print('reading slice beam', datafile)
                beam.read_HDF5_beam_file(datafile)
                # print('plotting slice beam', datafile)
                self.addsliceDataObject(beam, id=id, color=color)

    def changeSliceLengths(self):
        ''' change the time slice length for all beam objects '''
        for d in self.beams:
            self.changeSliceLength(d)
        # self.updateMultiAxisPlot()

    def changeSliceLength(self, name):
        ''' change the time slice length for a beam data object and update the plot '''
        beam = self.beams[name]
        beam.slices = self.slicePlotSliceWidthWidget.value()
        beam.bin_time()
        for n, param in enumerate(self.plotParams):
            if not param == 'next_row':
                label = param['label']
                exponent = np.floor(np.log10(np.abs(beam.slice_length)))
                x = 10**(12) * np.array((beam.slice_bins - np.mean(beam.slice_bins)))
                # self.multiPlotWidgets[label].setRange(xRange=[min(x),max(x)])
                y = getattr(beam, param['quantity'])
                if name in self.curves and label in self.curves[name]:
                    self.updateCurve(x, y, name, label)
                # self.curves[datafile][label].setData(x=x, y=y)

    def clear(self):
        self.beams = {}
        super(slicePlotWidget, self).clear()

def main():
    app = QApplication(sys.argv)
    pg.setConfigOptions(antialias=True)
    pg.setConfigOption('background', 'w')
    pg.setConfigOption('foreground', 'k')
    ex = slicePlotter()
    ex.show()
    ex.slicePlot.addsliceDataFiles([
    {'directory': 'OnlineModel_test_data/basefiles_4_250pC', 'filename': 'CLA-S02-APER-01.hdf5'},
    {'directory': 'OnlineModel_test_data/test_4', 'filename': ['CLA-L02-APER.hdf5','CLA-S04-APER-01.hdf5']}])
    ex.slicePlot.addsliceDataFile('OnlineModel_test_data/test_4/CLA-S03-APER.hdf5')
    # ex.multiPlot.removePlot('base+4')
    sys.exit(app.exec_())

if __name__ == '__main__':
   main()
