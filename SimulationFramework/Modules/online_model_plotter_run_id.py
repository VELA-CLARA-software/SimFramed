import sys, os, time, math, datetime, copy, re,  h5py
from collections import OrderedDict
# import glob
try:
    from PyQt4.QtCore import *
    from PyQt4.QtGui import *
except:
    from PyQt5.QtCore import *
    from PyQt5.QtGui import *
    from PyQt5.QtWidgets import *
import pyqtgraph as pg
from pyqtgraph.graphicsItems.LegendItem import ItemSample
# import numpy as np
sys.path.append(os.path.abspath(os.path.realpath(__file__)+'/../../../'))
# print (sys.path)
# import SimulationFramework.Modules.read_beam_file as raf
# import SimulationFramework.Modules.read_twiss_file as rtf
from SimulationFramework.Modules.online_model_twissPlot import twissPlotWidget
from SimulationFramework.Modules.online_model_slicePlot import slicePlotWidget
from SimulationFramework.Modules.online_model_beamPlot import beamPlotWidget

sys.path.append(os.path.realpath(__file__)+'/../../../../')

class online_model_plotter(QMainWindow):
    def __init__(self, screenpositions, parent = None, directory='.'):
        super(online_model_plotter, self).__init__(parent)
        self.resize(1200,900)

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
        exitAction.triggered.connect(self.close)
        fileMenu.addAction(exitAction)


        self.onlineModelPlotter = onlineModelPlotterWidget(screenpositions, directory=directory)
        self.setCentralWidget(self.onlineModelPlotter)

class onlineModelPlotterWidget(QWidget):

    colors = [QColor('#F5973A'),QColor('#A95AA1'),QColor('#85C059'),QColor('#0F2080'),QColor('#BDB8AD'), 'r', 'k', 'm', 'g']

    def __init__(self, screenpositions, parent = None, directory='.'):
        super(onlineModelPlotterWidget, self).__init__()
        self.screenpositions = screenpositions
        self.directory = directory
        self.layout = QGridLayout()
        self.setLayout(self.layout)

        self.tabWidget = QTabWidget()

        self.twissPlotWidget = twissPlotWidget()
        self.slicePlotWidget = slicePlotWidget(ymin=0)
        self.beamPlotWidget = beamPlotWidget()

        self.plotColor = 0
        self.run_id_prefixes = {}
        self.run_id_color = {}

        self.fileSelector = QComboBox()
        self.fileSelector.setMinimumWidth(200)
        self.fileSelector.currentIndexChanged.connect(self.loadScreen)
        self.plotScreenButton = QPushButton('Plot')
        self.plotScreenButton.clicked.connect(self.loadScreen)
        self.removeScreenButton = QPushButton('Remove')
        self.removeScreenButton.clicked.connect(self.clearBeamScreens)
        self.removeScreenButton.setEnabled(False)
        self.beamWidget = QGroupBox()
        self.beamLayout = QHBoxLayout()
        self.beamLayout.addWidget(self.fileSelector)
        # self.beamLayout.addWidget(self.screenSelector)
        # self.beamLayout.addWidget(self.plotScreenButton)
        # self.beamLayout.addWidget(self.removeScreenButton)
        self.beamWidget.setLayout(self.beamLayout)
        self.beamWidget.setMaximumWidth(800)
        # self.beamWidget.setVisible(True)

        self.sliceWidget = QGroupBox()
        self.sliceWidget.setVisible(False)
        self.sliceLayout = QHBoxLayout()
        self.sliceWidget.setLayout(self.sliceLayout)
        self.sliceLayout.addWidget(self.slicePlotWidget.slicePlotSliceWidthWidget)
        self.sliceWidget.setMaximumWidth(150)

        self.folderBeamWidget = QWidget()
        self.folderBeamLayout = QHBoxLayout()
        self.folderBeamLayout.setAlignment(Qt.AlignLeft);
        self.folderBeamWidget.setLayout(self.folderBeamLayout)
        # self.folderBeamLayout.addWidget(self.folderWidget)
        self.folderBeamLayout.addWidget(self.beamWidget)
        self.folderBeamLayout.addWidget(self.sliceWidget)
        self.plotType = 'Twiss'
        self.tabWidget.addTab(self.twissPlotWidget,'Twiss Plots')
        self.tabWidget.addTab(self.beamPlotWidget,'Beam Plots')
        self.tabWidget.addTab(self.slicePlotWidget,'Slice Beam Plots')

        self.listWidget = QListWidget()
        self.listWidget.setMaximumWidth(200)
        # self.listWidget.itemDoubleClicked.connect(self.removeScreen)
        self.listWidget.itemClicked.connect(self.curveClicked)

        self.layout.addWidget(self.folderBeamWidget,0,0,1,2)
        self.layout.addWidget(self.listWidget,1,0,1,1)
        self.layout.addWidget(self.tabWidget,1,1,1,1)

        self.twissDataCounter = 0
        self.changeDirectory(self.directory)
        self.shadowCurves = []
        self.connect_plot_signals()
        self.tabWidget.currentChanged.connect(self.changeTab)
        self.fileSelector.currentIndexChanged.connect(self.loadScreen)

    def connect_plot_signals(self):
        # When either subplot highlights a plot, connect it to the other plot and the listWidget
        self.twissPlotWidget.highlightCurveSignal.connect(self.subplotHighlighted)
        self.slicePlotWidget.highlightCurveSignal.connect(self.subplotHighlighted)
        self.beamPlotWidget.highlightCurveSignal.connect(self.subplotHighlighted)
        self.twissPlotWidget.unHighlightCurveSignal.connect(self.subplotUnHighlighted)
        self.slicePlotWidget.unHighlightCurveSignal.connect(self.subplotUnHighlighted)
        self.beamPlotWidget.unHighlightCurveSignal.connect(self.subplotUnHighlighted)

    def loadBeamDataFile(self, directory, beamFileName, color, id):
        if os.path.isfile(directory+'/'+beamFileName):
            self.beamPlotWidget.addbeamDataFile(directory, beamFileName, id=id, color=color)
            self.slicePlotWidget.addsliceDataFile(directory+'/'+beamFileName, id=id, color=color)

    def addRunIDToListWidget(self, run_id, prefixes, color=None):
        self.run_id_prefixes[run_id] = prefixes
        if color is None:
            color = self.colors[self.plotColor % len(self.colors)]
            self.plotColor += 1
        self.run_id_color[run_id] = color
        widgetItem = QListWidgetItem()
        widgetItem.setFont(QFont('Verdana', weight=QFont.Normal))
        widgetItem.setText(run_id)
        pixmap = QPixmap(16, 16)
        if not isinstance(color, QColor):
            color = pg.mkColor(color)
        self.run_id_color[run_id] = color
        pixmap.fill(color)
        icon = QIcon(pixmap)
        widgetItem.setIcon(icon)
        self.listWidget.addItem(widgetItem)
        self.loadTwiss(run_id)
        self.loadScreen()
        return color

    def removeRunIDFromListWidget(self, run_id):
        items = self.listWidget.findItems(run_id, Qt.MatchExactly)
        if len(items) > 0:
            item = items[0]
            row = self.listWidget.row(item)
            self.listWidget.takeItem(row)
        del self.run_id_prefixes[run_id]
        self.twissPlotWidget.removeCurve(run_id)
        self.slicePlotWidget.removeCurve(run_id)
        self.beamPlotWidget.removePlot(run_id)

    def updatePlotButtons(self):
        if len(self.screenpositions) > 0 and not self.fileSelector.currentText() == '':
            self.beamFileName = str(self.fileSelector.currentData())+'.hdf5'
            if self.directory+'/'+self.beamFileName in self.slicePlotWidget.curves:
                self.plotScreenButton.setText('Update')
                self.removeScreenButton.setEnabled(True)
            else:
                self.plotScreenButton.setText('Plot')
                self.removeScreenButton.setEnabled(False)

    def changeTab(self, i):
        if self.tabWidget.tabText(i) == 'Beam Plots':
            self.plotType = 'Beam'
            # self.beamWidget.setVisible(True)
            self.sliceWidget.setVisible(False)
        elif self.tabWidget.tabText(i) == 'Slice Beam Plots':
            self.plotType = 'Slice'
            # self.beamWidget.setVisible(True)
            self.sliceWidget.setVisible(True)
        else:
            self.plotType = 'Twiss'
            # self.beamWidget.setVisible(False)
            self.sliceWidget.setVisible(False)

    def changeDirectory(self, directory=None, id=None):
        # pass
        # self.directory = directory
        self.currentFileText = self.fileSelector.currentText()
        self.updateFileCombo()
        # self.loadTwissDataFile()

    def updateFileCombo(self):
        self.fileSelector.clear()
        i = -1
        self.allscreens = []
        for l, f in self.screenpositions.items():
            for k, v in f.items():
                n = k
                p = v['position']
                self.allscreens.append([n, p, l])
        sortedscreennames = sorted(self.allscreens, key=lambda x: float(x[1]))
        selected = False
        for n,p,l in sortedscreennames:
            self.fileSelector.addItem(n.ljust(20,' ') + '('+str(p)+'m)',n)
            i += 1
            if n == self.currentFileText:
                self.fileSelector.setCurrentIndex(i)
                selected = True
        if not selected:
            self.fileSelector.setCurrentIndex(3)

    def loadTwiss(self, id):
        prefixes = self.run_id_prefixes[id]
        twissList = []
        for s, d in prefixes.items():
            twissList.append({'directory': 'test/'+d, 'sections': [s]})
        print('twissList = ', twissList)
        self.twissPlotWidget.addTwissDirectory(twissList, id=id, color=self.run_id_color[id])

    def loadScreen(self):
        self.clearBeamScreens()
        if len(self.screenpositions) > 0 and not self.fileSelector.currentText() == '':
            beamfilename = str(self.fileSelector.currentData())+'.hdf5'
            screens, positions, lattices = list(zip(*self.allscreens))
            screen_idx = screens.index(self.fileSelector.currentData())
            lattice = lattices[screen_idx]
            for run, prefixes in self.run_id_prefixes.items():
                directory = 'test/' + prefixes[lattice]
                color = self.run_id_color[run]
                self.loadBeamDataFile(directory, beamfilename, color, run)

    def curveClicked(self, item):
        name = item.text()
        if not name in self.shadowCurves:
            self.highlightPlot(item)
        else:
            self.unHighlightPlot(item)

    def get_run_id_directory(self, run_id):
        beamfilename = str(self.fileSelector.currentData())+'.hdf5'
        screens, positions, lattices = list(zip(*self.allscreens))
        screen_idx = screens.index(self.fileSelector.currentData())
        lattice = lattices[screen_idx]
        prefix = self.run_id_prefixes[run_id]
        directory = 'test/' + prefix[lattice]
        return directory+'/'+run_id

    def highlightPlot(self, item):
        name = item.text()
        if not name in self.shadowCurves:
            self.shadowCurves.append(name)
            item.setFont(QFont('Verdana', weight=QFont.Bold))
        self.twissPlotWidget.highlightPlot(name)
        self.slicePlotWidget.highlightPlot(name)
        self.beamPlotWidget.highlightPlot(name)

    def unHighlightPlot(self, item):
        name = item.text()
        if name in self.shadowCurves:
            self.shadowCurves.remove(name)
            item.setFont(QFont('Verdana', weight=QFont.Normal))
        self.twissPlotWidget.unHighlightPlot(name)
        self.slicePlotWidget.unHighlightPlot(name)
        self.beamPlotWidget.unHighlightPlot(name)

    def subplotHighlighted(self, name):
        subname = name.split('/')[-1]
        items = self.listWidget.findItems(subname, Qt.MatchExactly)
        if len(items) > 0:
            self.highlightPlot(items[0])

    def subplotUnHighlighted(self, name):
        subname = name.split('/')[-1]
        items = self.listWidget.findItems(subname, Qt.MatchExactly)
        if len(items) > 0:
            self.unHighlightPlot(items[0])

    def clearBeamScreens(self):
        self.beamPlotWidget.clear()
        self.slicePlotWidget.clear()

def main():
    # global app
    import argparse
    parser = argparse.ArgumentParser(description='Analyse Online Model Folder')
    parser.add_argument('-d', '--directory', default='.', type=str)
    args = parser.parse_args()
    import ruamel.yaml as yaml
    yaml.add_representer(OrderedDict, yaml.representer.SafeRepresenter.represent_dict)
    with open(r'C:\Users\jkj62\Documents\GitHub\ASTRA_COMPARISONRunner-HMCC\screen_positions.yaml', 'r') as stream:
        yaml_parameter_dict = yaml.safe_load(stream)
    # print(yaml_parameter_dict)
    pg.setConfigOptions(antialias=True)
    pg.setConfigOption('background', 'w')
    pg.setConfigOption('foreground', 'k')
    app = QApplication(sys.argv)
    ex = online_model_plotter(yaml_parameter_dict,directory=args.directory)
    ex.show()
    sys.exit(app.exec_())

if __name__ == '__main__':
   main()
