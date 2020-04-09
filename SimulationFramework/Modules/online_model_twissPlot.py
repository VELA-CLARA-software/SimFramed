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
    ''' QWidget containing various Twiss plots '''
    # Styles for the plot lines
    colors = [QColor('#F5973A'),QColor('#A95AA1'),QColor('#85C059'),QColor('#0F2080'),QColor('#BDB8AD'), 'r']
    styles = [Qt.SolidLine, Qt.DashLine, Qt.DotLine, Qt.DashDotLine, Qt.DashDotDotLine]

    # Layout oder for the individual Tiwss plot items
    twissParams = [
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
        ''' twissPlotWidget - main pyQtGraph display widgets '''
        self.twissPlotView = pg.GraphicsView(useOpenGL=True)
        self.twissPlotWidget = pg.GraphicsLayout()
        self.twissPlotView.setCentralItem(self.twissPlotWidget)
        ''' twissPlotWidgets - holds the base plotWidgets for each Twiss plot '''
        self.twissPlotWidgets = {}
        ''' curves - a dictionary containing {directory: [individual plotItems for each twiss plot,]} '''
        self.curves = {}
        for n, param in enumerate(self.twissParams):
            if param == 'next_row':
                self.twissPlotWidget.nextRow()
            else:
                ''' p - the relevant plotWidget for each param in twissParams '''
                p = self.twissPlotWidget.addPlot(title=param['label'])
                ''' this just links the horizontal axis for all plots.
                    The first plot is selected to be the master axis '''
                if n == 0:
                    self.linkAxis = p.vb
                else:
                    p.setXLink(self.linkAxis)
                p.showGrid(x=True, y=True)
                p.setLabel('left', text=param['name'], units=param['units'])
                ''' set lower plot limit at 0 for all plots '''
                p.vb.setLimits(xMin=0,yMin=0)
                ''' set the vertical viewRange according to the twissParams '''
                p.vb.setYRange(*param['range'])
                ''' record the plotWidget in the dictionary indexed by Twiss item '''
                self.twissPlotWidgets[param['label']] = p

        self.layout = QVBoxLayout()
        self.setLayout(self.layout)

        self.layout.addWidget(self.twissPlotView)

        ''' used for style cycling '''
        self.plotColor = 0
        self.shadowCurves = []

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
        twiss = self.loadDataFile(**directory[0], reset=True)
        ''' loads any other directories '''
        for d in directory[1:]:
            twiss = self.loadDataFile(**d, reset=False, twiss=twiss)
        ''' assignes a reference name if none is given '''
        name = directory[-1]['directory'] if name is None else name
        self.addtwissDataObject(twiss, name)

    def addtwissDataObject(self, dataobject, indexname):
        '''
            Take a twiss data object and add a plot line to all of the relevant Twiss plots

            Keyword arguments:
            dataobject -- twiss data object for plotting
            indexname -- dictionary key name

        '''
        twiss = dataobject
        ''' select a color and style '''
        color = self.colors[self.plotColor % len(self.colors)]
        pen = pg.mkPen(color, width=3, style=self.styles[int(self.plotColor / len(self.styles))])
        ''' iterate the color index '''
        self.plotColor += 1

        self.curves[indexname] = {}
        for param in self.twissParams:
            if not param == 'next_row':
                label = param['label']
                if len(twiss[param['quantity']]) > 0: # confirm the data is there!
                    ''' load the data in z (ASTRA does not do s-coordinates) and then ensure it is sorted correctly '''
                    x = twiss['z']
                    y = twiss[param['quantity']]
                    xy = np.transpose(np.array([x,y]))
                    x, y = np.transpose(xy[np.argsort(xy[:,0])])
                    ''' create a plotItem on a Twiss plot and save to the curves dictionary '''
                    self.curves[indexname][label] = self.twissPlotWidgets[label].plot(x=x, y=y, pen=pen)
                    self.curves[indexname][label].curve.setClickable(True)
                    self.curves[indexname][label].sigClicked.connect(lambda: self.highlightPlot(indexname))

    def removePlot(self, name):
        ''' finds all Twiss plots based on a key name, and removes them '''
        if not isinstance(name, (list, tuple)):
            name = [name]
        for n in name:
            if n in self.curves:
                for param in self.twissParams:
                    if not param == 'next_row':
                        ''' Remove the plotItem from the relevant plotWidget '''
                        self.twissPlotWidgets[param['label']].removeItem(self.curves[n][param['label']])

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

    def highlightPlot(self, name):
        ''' highlights a particular plot '''
        # print('highligher clicked! = ', name)
        if not isinstance(name, (list, tuple)):
            name = [name]
        for n in name:
            self.addShadowPen(n)
            for n in self.curves.keys():
                if n in self.shadowCurves or not len(self.shadowCurves) > 0:
                    self.setPenAlpha(n, 255, 3)
                else:
                    self.setPenAlpha(n, 25, 3)

    def addShadowPen(self, name):
        for param in self.twissParams:
            if not param == 'next_row':
                label = param['label']
                curve = self.curves[name][label]
                if curve.opts['shadowPen'] is None:
                    self.shadowCurves.append(name)
                    pen = curve.opts['pen']
                    shadowpencolor = pen.color()
                    shadowpencolor.setAlpha(100)
                    shadowpen = pg.mkPen(color=shadowpencolor, width=(pen.width()+3))
                    curve.setShadowPen(shadowpen)
                else:
                    self.shadowCurves.remove(name)
                    curve.setShadowPen(None)
                    curve.opts['shadowPen'] = None

    def setPenAlpha(self, name, alpha=255, width=3):
        for param in self.twissParams:
            if not param == 'next_row':
                label = param['label']
                curve = self.curves[name][label]
                pen = curve.opts['pen']
                pencolor = pen.color()
                pencolor.setAlpha(alpha)
                pen = pg.mkPen(color=pencolor, width=width, style=pen.style())
                curve.setPen(pen)


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
    ex.twissPlot.addTwissDirectory([{'directory': 'OnlineModel_test_data/basefiles_4_250pC', 'sections': ['injector400']}, {'directory': 'OnlineModel_test_data/test_4', 'sections': 'All'}], name='base+4')
    ex.twissPlot.addTwissDirectory([{'directory': 'OnlineModel_test_data/basefiles_4_250pC', 'sections': ['injector400']}, {'directory': 'OnlineModel_test_data/test_2', 'sections': 'All'}], name='base+2')
    # ex.twissPlot.removePlot('base+4')
    sys.exit(app.exec_())

if __name__ == '__main__':
   main()
