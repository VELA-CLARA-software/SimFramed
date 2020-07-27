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
from SimulationFramework.Modules.multiPlot import multiPlotWidget
sys.path.append(os.path.realpath(__file__)+'/../../../../')

class mainWindow(QMainWindow):
    def __init__(self):
        super(mainWindow, self).__init__()
        self.resize(1800,900)
        self.centralWidget = QWidget()
        self.layout = QVBoxLayout()
        self.centralWidget.setLayout(self.layout)

        self.tab = QTabWidget()
        self.twissPlotWidget = twissPlotWidget()

        self.layout.addWidget(self.twissPlotWidget)

        self.setCentralWidget(self.centralWidget)

        self.setWindowTitle("ASTRA Data Plotter")
        menubar = self.menuBar()
        fileMenu = menubar.addMenu('&File')

        exitAction = QAction('&Exit', self)
        exitAction.setShortcut('Ctrl+Q')
        exitAction.setStatusTip('Exit application')
        exitAction.triggered.connect(self.close)
        fileMenu.addAction(exitAction)

class twissPlotWidget(multiPlotWidget):
    # Layout oder for the individual Tiwss plot items
    plotParams = [
                   {'label': 'Horizontal Beam Size', 'name': '&sigma;<sub>x</sub>', 'quantity': 'sigma_x', 'range': [0,1e-3], 'units': 'm'},
                   {'label': 'Vertical Beam Size', 'name': '&sigma;<sub>y</sub>', 'quantity': 'sigma_y', 'range': [0,1e-3], 'units': 'm'},
                   {'label': 'Momentum', 'name': 'cp', 'quantity': 'cp_eV', 'range': [0,250e6], 'units': 'eV/c'},
                   'next_row',
                   {'label': 'Momentum Spread', 'name': '&sigma;<sub>cp</sub>', 'quantity': 'sigma_cp_eV', 'range': [0,5e6], 'units': 'eV/c'},
                   {'label': 'Bunch Length', 'name': '&sigma;<sub>z</sub>', 'quantity': 'sigma_z', 'range': [0,0.6e-3], 'units': 'm'},
                   {'label': 'Horizontal Emittance (normalised)', 'name': '&epsilon;<sub>n,x</sub>', 'quantity': 'enx', 'range': [0.,1.5e-6], 'units': 'm-rad'},
                   'next_row',
                   {'label': 'Vertical Emittance (normalised)', 'name': '&epsilon;<sub>n,y</sub>', 'quantity': 'eny', 'range':  [0.,1.5e-6], 'units': 'm-rad'},
                   {'label': 'Horizontal Beta Function', 'name': '&beta;<sub>x</sub>', 'quantity': 'beta_x', 'range': [0,200], 'units': 'm'},
                   {'label': 'Vertical Beta Function', 'name': '&beta;<sub>y</sub>', 'quantity': 'beta_y', 'range': [0,200], 'units': 'm'},
                  ]

    def __init__(self, **kwargs):
        super(twissPlotWidget, self).__init__(xmin=0, ymin=0, **kwargs)

    def addTwissDirectory(self, directory, name=None):
        '''
            Read the data files in a directory and add a plotItem to the relevant curves

            Keyword arguments:
            directory -- dictionary containing directory definitions:
                [
                    {'directory': <dir location>,           'sections': [<list of lattice sections>]},
                    {'directory': <another dir location>,   'sections': 'All'},
                    ...
                ]
                The directories are scanned for ASTRA or Elegant twiss files and imported as one set.
                The data is then sorted in z. No analysis is done for duplicate entries.
            name -- dictionary key name (default: last directory name)
        '''
        ''' load the data files into the twiss dictionary '''
        if not isinstance(directory, (list, tuple)):
            directory = [directory]
        ''' loads the first (and only?) param in the list of directories '''
        twiss = self.loadDataFile( reset=True, **(directory[0]))
        ''' loads any other directories '''
        for d in directory[1:]:
            twiss = self.loadDataFile(reset=False, twiss=twiss, **d)
        ''' assignes a reference name if none is given '''
        name = directory[-1]['directory'] if name is None else name
        self.addtwissDataObject(twiss, name)

    def addtwissDataObject(self, dataobject, name):
        '''
            Take a twiss data object and add a plot line to all of the relevant Twiss plots

            Keyword arguments:
            dataobject -- twiss data object for plotting
            name -- dictionary key name

        '''
        twiss = dataobject
        ''' select a color and style '''
        color = self.colors[self.plotColor % len(self.colors)]
        style = self.styles[int(self.plotColor % len(self.styles))]
        pen = pg.mkPen(color, width=3, style=style)
        ''' iterate the color index '''
        if name not in self.curves:
            self.curves[name] = {}
            self.plotColor += 1
        for param in self.plotParams:
            if not param == 'next_row':
                label = param['label']
                if len(twiss[param['quantity']]) > 0: # confirm the data is there!
                    ''' load the data in z (ASTRA does not do s-coordinates) and then ensure it is sorted correctly '''
                    x = twiss['z']
                    y = twiss[param['quantity']]
                    xy = np.transpose(np.array([x,y]))
                    x, y = np.transpose(xy[np.argsort(xy[:,0])])
                    ''' create a plotItem on a Twiss plot and save to the curves dictionary '''
                    if name in self.curves and label in self.curves[name]:
                        print('Updating curve: ', name, label)
                        self.updateCurve(x, y, name, label)
                    else:
                        print('ADDING curve: ', name, label)
                        self.addCurve(x, y, name, label, pen)
        if not isinstance(color, QColor):
            color = pg.mkColor(color)
        return color, style

    def loadDataFile(self, directory, sections=None, reset=True, twiss=None):
        ''' loads ASTRA and Elegant data files from a directory and returns a twiss object'''
        if sections is None or (isinstance(sections, str) and sections.lower() == 'all'):
            astrafiles = sorted(glob.glob(directory+"/*Xemit*"))
            elegantfiles = sorted(glob.glob(directory+"/*.flr"))
        else:
            astrafiles = []
            elegantfiles = []
            for s in sections:
                astrafiles += sorted(glob.glob(directory+"/"+s+"*Xemit*"))
                elegantfiles += sorted(glob.glob(directory+"/"+s+"*.flr"))
        if twiss is None: # If it doesn't exist need to instantiate a twiss obkject
            twiss = rtf.twiss()
        twiss.read_astra_emit_files(astrafiles, reset=reset)
        reset = False if len(astrafiles) > 0 else reset # if we have alreay found some ASTRA files, we need to set this to false to append new data, otherwise check input value
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
    ex.twissPlotWidget.addTwissDirectory([{'directory': 'OnlineModel_test_data/basefiles_4_250pC', 'sections': ['injector400']}, {'directory': 'OnlineModel_test_data/test_4', 'sections': 'All'}], name='base+4')
    ex.twissPlotWidget.addTwissDirectory([{'directory': 'OnlineModel_test_data/basefiles_4_250pC', 'sections': ['injector400']}, {'directory': 'OnlineModel_test_data/test_2', 'sections': 'All'}], name='base+2')
    # ex.multiPlot.removePlot('base+4')
    sys.exit(app.exec_())

if __name__ == '__main__':
   main()
