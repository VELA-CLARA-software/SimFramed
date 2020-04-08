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
    colors = [QColor('#F5973A'),QColor('#A95AA1'),QColor('#85C059'),QColor('#0F2080'),QColor('#BDB8AD'), 'r', 'k', 'm', 'g']
    styles = [Qt.SolidLine, Qt.DashLine, Qt.DotLine, Qt.DashDotLine, Qt.DashDotDotLine]

    twissParams = [
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
        ''' These are for reading data files from ASTRA and Elegant '''
        self.twissObjects = {}
        # self.twiss = rtf.twiss()

        self.twissPlotWidget = QWidget()
        self.twissPlotLayout = QVBoxLayout()
        self.twissPlotWidget.setLayout(self.twissPlotLayout)
        self.twissPlotWidgetGraphicsLayout = pg.GraphicsLayoutWidget()
        self.twissPlotCheckbox = {}
        self.viewboxes = {}
        self.curves = {}
        self.twissaxis = {}
        self.twissPlotCheckboxWidget = QWidget()
        self.twissPlotCheckboxLayout = QVBoxLayout()
        self.twissPlotCheckboxWidget.setLayout(self.twissPlotCheckboxLayout)
        self.twissPlot = self.twissPlotWidgetGraphicsLayout.addPlot(title='Twiss',row=0,col=50)
        self.twissPlot.showAxis('left', False)
        self.twissPlot.showGrid(x=True, y=True)
        self.twissPlot.enableAutoRange(x=True, y=True)
        self.twissPlot.vb.setLimits(xMin=0,yMin=0)
        i = 0
        for param in self.twissParams:
            axis = pg.AxisItem("left")
            labelStyle = {'color': '#'+pg.colorStr(pg.mkColor(self.colors[i]))[0:-2]}
            axis.setLabel(text=param['name'], units=param['units'], **labelStyle)
            i += 1
            viewbox = pg.ViewBox()
            self.viewboxes[param['label']] = viewbox
            axis.linkToView(viewbox)
            viewbox.setXLink(self.twissPlot.vb)
            self.twissaxis[param['label']] = [axis, viewbox]
            self.curves[param['label']] = {}
            col = i
            self.twissPlotWidgetGraphicsLayout.ci.addItem(axis, row = 0, col = col,  rowspan=1, colspan=1)
            self.twissPlotWidgetGraphicsLayout.ci.addItem(viewbox, row=0, col=50)
            viewbox.setLimits(yMin=0)
            self.twissPlotCheckbox[param['label']] = QCheckBox(param['label'])
            self.twissPlotCheckbox[param['label']].setChecked(True)
            self.twissPlotCheckboxLayout.addWidget(self.twissPlotCheckbox[param['label']])
            self.twissPlotCheckbox[param['label']].stateChanged.connect(self.updateTwissPlot)

        self.twissPlotAxisWidget = QWidget()
        self.twissPlotAxisLayout = QHBoxLayout()
        self.twissPlotAxisWidget.setLayout(self.twissPlotAxisLayout)
        self.twissPlotAxisLayout.addWidget(self.twissPlotCheckboxWidget)
        self.twissPlotLayout.addWidget(self.twissPlotAxisWidget)
        self.twissPlotLayout.addWidget(self.twissPlotWidgetGraphicsLayout)

        self.layout = QVBoxLayout()
        self.setLayout(self.layout)

        self.layout.addWidget(self.twissPlotWidget)

        ''' used for style cycling '''
        self.plotColor = 0

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

    def addtwissDataObject(self, twissobject, datafile):
        ''' addtwissDirectory - read the data files in a directory and add a plotItem to the relevant twissPlotItems '''
        ''' load the data files into the twiss dictionary '''
        if str(type(twissobject)) == "<class 'SimulationFramework.Modules.read_twiss_file.twiss'>":
            for n, param in enumerate(self.twissParams):
                color = self.colors[n]
                pen = pg.mkPen(color=color, style=self.styles[self.plotColor % len(self.styles)])
                self.curves[param['label']][datafile] = pg.PlotDataItem()
                self.viewboxes[param['label']].addItem(self.curves[param['label']][datafile])
                x = twissobject['z']
                y = twissobject[param['quantity']]
                xy = np.transpose(np.array([x,y]))
                x, y = np.transpose(xy[np.argsort(xy[:,0])])
                self.curves[param['label']][datafile].setData(x=x, y=y, pen=pen)
            self.plotColor += 1
        self.updateTwissPlot()

    def updateTwissPlot(self):
        for n, param in enumerate(self.twissParams):
            if self.twissPlotCheckbox[param['label']].isChecked():
                for d in self.curves[param['label']]:
                    self.curves[param['label']][d].setVisible(True)
                self.twissaxis[param['label']][0].setVisible(True)
            else:
                for d in self.curves[param['label']]:
                    self.curves[param['label']][d].setVisible(False)
                self.twissaxis[param['label']][0].setVisible(False)
            self.twissaxis[param['label']][1].autoRange()
            currentrange = self.twissaxis[param['label']][1].viewRange()
            self.twissaxis[param['label']][1].setXRange(0, currentrange[0][1])
            self.twissaxis[param['label']][1].setYRange(0, currentrange[1][1])

    def removePlot(self, directory):
        ''' finds all twiss plots based on a directory name, and removes them '''
        if not isinstance(directory, (list, tuple)):
            directory = [directory]
        indexname = directory[-1]
        if directory in self.twissPlotItems:
            for entry in self.twissplotLayout:
                if entry == 'next_row':
                    pass
                else:
                    self.twissPlotWidgets[entry['label']].removeItem(self.twissPlotItems[indexname][entry['label']])

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
