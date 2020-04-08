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

class twissPlotWidget(QWidget):
    # Styles for the plot lines
    colors = [QColor('#F5973A'),QColor('#A95AA1'),QColor('#85C059'),QColor('#0F2080'),QColor('#BDB8AD'), 'r']
    styles = [Qt.SolidLine, Qt.DashLine, Qt.DotLine, Qt.DashDotLine, Qt.DashDotDotLine]

    # Layout oder for the individual Tiwss plot items
    twissplotLayout = [
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
        super(twissPlotWidget, self).__init__(**kwargs)
        ''' These are for reading data files from ASTRA and Elegant '''
        self.beam = raf.beam()
        self.twiss = rtf.twiss()

        ''' twissPlotWidget - main pyQtGraph display widgets '''
        self.twissPlotView = pg.GraphicsView(useOpenGL=True)
        self.twissPlotWidget = pg.GraphicsLayout()
        self.twissPlotView.setCentralItem(self.twissPlotWidget)
        ''' twissPlotWidgets - holds the base plotWidgets for each Twiss plot '''
        self.twissPlotWidgets = {}
        ''' twissPlotItems - a dictionary containing {directory: [individual plotItems for each twiss plot,]} '''
        self.twissPlotItems = {}
        for n, entry in enumerate(self.twissplotLayout):
            if entry == 'next_row':
                self.twissPlotWidget.nextRow()
            else:
                ''' p - the relevant plotWidget for each entry in twissplotLayout '''
                p = self.twissPlotWidget.addPlot(title=entry['label'])
                ''' this just links the horizontal axis for all plots.
                    The first plot is selected to be the master axis '''
                if n == 0:
                    self.linkAxis = p.vb
                else:
                    p.setXLink(self.linkAxis)
                p.showGrid(x=True, y=True)
                p.setLabel('left', text=entry['name'], units=entry['units'])
                ''' vb = current plots viewbox '''
                vb = p.vb
                ''' set lower plot limit at 0 for all plots '''
                p.vb.setLimits(xMin=0,yMin=0)
                ''' set the vertical viewRange according to the twissplotLayout '''
                vb.setYRange(*entry['range'])
                ''' record the plotWidget in the dictionary indexed by Twiss item '''
                self.twissPlotWidgets[entry['label']] = p

        self.layout = QVBoxLayout()
        self.setLayout(self.layout)

        self.layout.addWidget(self.twissPlotView)

        ''' used for style cycling '''
        self.plotColor = 0

    def addTwissDirectory(self, directory):
        ''' addTwissDirectory - read the data files in a directory and add a plotItem to the relevant twissPlotItems '''
        ''' select a color and style '''
        color = self.colors[self.plotColor % len(self.colors)]
        pen = pg.mkPen(color, width=3, style=self.styles[int(self.plotColor / len(self.styles))])
        self.plotColor += 1
        ''' load the data files into the twiss dictionary '''
        datadict=[]
        if not isinstance(directory, (list, tuple)):
            directory = [directory]
        self.loadDataFile(**directory[0], reset=True)
        for d in directory[1:]:
            self.loadDataFile(**d, reset=False)
        indexname = directory[-1]['directory']
        self.twissPlotItems[indexname] = {}
        for entry in self.twissplotLayout:
            if entry == 'next_row':
                pass
            else:
                if len(self.twiss[entry['quantity']]) > 0:
                    ''' load the data in z (ASTRA does not s-coordinates) and then ensure it is sorted correctly '''
                    x = self.twiss['z']
                    y = self.twiss[entry['quantity']]
                    xy = np.transpose(np.array([x,y]))
                    x, y = np.transpose(xy[np.argsort(xy[:,0])])
                    ''' create a plotItem on a Twiss plot and save to the twissPlotItems dictionary '''
                    self.twissPlotItems[indexname][entry['label']] = self.twissPlotWidgets[entry['label']].plot(x=x, y=y, pen=pen)

    def removePlot(self, directory):
        ''' finds all Twiss plots based on a directory name, and removes them '''
        if not isinstance(directory, (list, tuple)):
            directory = [directory]
        indexname = directory[-1]
        if directory in self.twissPlotItems:
            for entry in self.twissplotLayout:
                if entry == 'next_row':
                    pass
                else:
                    self.twissPlotWidgets[entry['label']].removeItem(self.twissPlotItems[indexname][entry['label']])

    def loadDataFile(self, directory, sections=None, reset=True):
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
        self.twiss.read_astra_emit_files(astrafiles, reset=reset)
        reset = False if len(astrafiles) > 0 else reset
        ''' reset=False stops the previously loaded data from being overwritten'''
        self.twiss.read_elegant_twiss_files(elegantfiles, reset=reset)

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
