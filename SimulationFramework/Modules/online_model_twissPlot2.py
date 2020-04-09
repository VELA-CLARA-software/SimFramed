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
from SimulationFramework.Modules.online_model_multiaxis_Plot import multiaxisPlotWidget
sys.path.append(os.path.realpath(__file__)+'/../../../../')

class mainWindow(QMainWindow):
    def __init__(self):
        super(mainWindow, self).__init__()
        self.resize(1800,900)
        self.centralWidget = QWidget()
        self.layout = QVBoxLayout()
        self.centralWidget.setLayout(self.layout)

        self.tab = QTabWidget()
        self.twissPlot = twissPlotWidget()

        self.layout.addWidget(self.twissPlot)

        self.setCentralWidget(self.centralWidget)

        self.setWindowTitle("ASTRA Data Plotter")
        menubar = self.menuBar()
        fileMenu = menubar.addMenu('&File')

        exitAction = QAction('&Exit', self)
        exitAction.setShortcut('Ctrl+Q')
        exitAction.setStatusTip('Exit application')
        exitAction.triggered.connect(self.close)
        fileMenu.addAction(exitAction)

class twissPlotWidget(multiaxisPlotWidget):

    plotParams = [
       {'label': 'Horizontal Beam Size', 'name': '&sigma;<sub>x</sub>', 'quantity': 'sigma_x', 'range': [0,1e-3], 'units': 'm'},
       {'label': 'Vertical Beam Size', 'name': '&sigma;<sub>y</sub>', 'quantity': 'sigma_y', 'range': [0,1e-3], 'units': 'm'},
       {'label': 'Momentum', 'name': 'cp', 'quantity': 'cp_eV', 'range': [0,250e6], 'units': 'eV/c'},
       {'label': 'Momentum Spread', 'name': '&sigma;<sub>cp</sub>', 'quantity': 'sigma_cp_eV', 'range': [0,5e6], 'units': 'eV/c'},
       {'label': 'Bunch Length', 'name': '&sigma;<sub>z</sub>', 'quantity': 'sigma_z', 'range': [0,0.6e-3], 'units': 'm'},
       {'label': 'Horizontal Emittance (normalised)', 'name': '&epsilon;<sub>n,x</sub>', 'quantity': 'enx', 'range': [0.,1.5e-6], 'units': 'm-rad'},
       {'label': 'Vertical Emittance (normalised)', 'name': '&epsilon;<sub>n,y</sub>', 'quantity': 'eny', 'range':  [0.,1.5e-6], 'units': 'm-rad'},
       {'label': 'Horizontal Beta Function', 'name': '&beta;<sub>x</sub>', 'quantity': 'beta_x', 'range': [0,200], 'units': 'm'},
       {'label': 'Vertical Beta Function', 'name': '&beta;<sub>y</sub>', 'quantity': 'beta_y', 'range': [0,200], 'units': 'm'},
    ]


    def __init__(self, **kwargs):
        super(twissPlotWidget, self).__init__(**kwargs)

    def addTwissDirectory(self, directory):
        ''' addTwissDirectory - read the data files in a directory and add a plotItem to the relevant twissPlotItems '''
        ''' load the data files into the twiss dictionary '''
        datadict=[]
        if not isinstance(directory, (list, tuple)):
            directory = [directory]
        twiss = self.loadDataFile(**directory[0], reset=True)
        for d in directory[1:]:
            twiss = self.loadDataFile(**d, reset=False, twiss=twiss)
        indexname = directory[-1]['directory']
        self.addtwissDataObject(twiss, directory[-1]['directory'])

    def addtwissDataObject(self, twissobject, name):
        ''' addtwissDirectory - read the data files in a directory and add a plotItem to the relevant twissPlotItems '''
        ''' load the data files into the twiss dictionary '''
        if str(type(twissobject)) == "<class 'SimulationFramework.Modules.read_twiss_file.twiss'>":
            self.curves[name] = {}
            for n, param in enumerate(self.plotParams):
                label = param['label']
                color = self.colors[n]
                pen = pg.mkPen(color=color, style=self.styles[self.plotColor % len(self.styles)])
                self.curves[name][label] = pg.PlotDataItem()
                self.curves[name][label].curve.setClickable(True)
                self.curves[name][label].sigClicked.connect(lambda: self.highlightPlot(name))
                self.viewboxes[label].addItem(self.curves[name][label])
                x = twissobject['z']
                y = twissobject[param['quantity']]
                xy = np.transpose(np.array([x,y]))
                x, y = np.transpose(xy[np.argsort(xy[:,0])])
                self.curves[name][label].setData(x=x, y=y, pen=pen)
            self.plotColor += 1
        self.updateMultiAxisPlot()

    def loadDataFile(self, directory, sections=None, reset=True, twiss=None):
        # print('loading directory = ', directory)
        ''' loads ASTRA and Elegant data files from a directory '''
        if sections is None or (isinstance(sections, str) and sections.lower() == 'all'):
            astrafiles = sorted(glob.glob(directory+"/*Xemit*"))
            elegantfiles = sorted(glob.glob(directory+"/*.flr"))
        else:
            astrafiles = []
            elegantfiles = []
            for s in sections:
                print(s)
                astrafiles += sorted(glob.glob(directory+"/"+s+"*Xemit*"))
                elegantfiles += sorted(glob.glob(directory+"/"+s+"*.flr"))
        if twiss is None:
            twiss = rtf.twiss()
        twiss.read_astra_emit_files(astrafiles, reset=reset)
        reset = False if len(astrafiles) > 0 else reset
        ''' reset=False stops the previously loaded data from being overwritten'''
        twiss.read_elegant_twiss_files(elegantfiles, reset=reset)
        return twiss

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
    ex.twissPlot.addTwissDirectory([{'directory': 'OnlineModel_test_data/basefiles_4_250pC', 'sections': ['injector400']}, {'directory': 'OnlineModel_test_data/test_4', 'sections': 'All'}])
    ex.twissPlot.addTwissDirectory([{'directory': 'OnlineModel_test_data/basefiles_4_250pC', 'sections': ['injector400']}, {'directory': 'OnlineModel_test_data/test_2', 'sections': 'All'}])
    sys.exit(app.exec_())

if __name__ == '__main__':
   main()
