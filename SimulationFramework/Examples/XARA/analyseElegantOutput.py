import sys
import time
import os
import csv
sys.path.append('../../../../S2E_FEL_dev_gen/')
from ocelot.adaptors.genesis import read_out_file
from ocelot.gui.genesis_plot import fwhm3
from ocelot.common.math_op import *
from scipy.stats import kde
import numpy as np
import argparse
sys.path.append(os.path.abspath(__file__+'/../../../../'))
import SimulationFramework.Modules.read_beam_file as rbf
import SimulationFramework.Modules.read_twiss_file as rtf
sys.path.append(os.path.abspath(__file__+'/../../../../../'))
import Software.Procedures.qt as qt
import pyqtgraph as pg
from pyqtgraph.widgets.MatplotlibWidget import MatplotlibWidget

parser = argparse.ArgumentParser(description='Analyse genesis output file and return relevant parameters')
parser.add_argument('-d', '--directory', default='.')
parser.add_argument('-f', '--file', default='run.0.gout')
parser.add_argument('-p', '--plot', default=True)

import matplotlib.pyplot as plt
import matplotlib.lines as mlines

def analyse_image(mw, dir=None):
    dir = str(dir)
    if dir is not None:
        if os.path.isfile(dir + '/' + 'CLA-S07-APER-01.hdf5'):
            args = parser.parse_args()
            beam = rbf.beam()
            twiss = rtf.twiss()
            beam.read_HDF5_beam_file(dir+'/CLA-S07-APER-01.hdf5')
            twiss.read_elegant_twiss_files(dir+'/CLARAX.twi' )
            beam.slice_length = 0.01e-12

            t = 1e12*(beam.t-np.mean(beam.t))
            t_grid = np.linspace(min(t), max(t), 2**8)
            # t_grid = np.arange(min(t), max(t), 0.01)
            peakIPDF = beam.PDFI(t, t_grid, bandwidth=beam.rms(t)/(2**3))*250
            peakICDF = beam.CDF(t, t_grid, bandwidth=beam.rms(t)/(2**3))
            peakIFWHM, indexes = beam.FWHM(t_grid, peakIPDF, frac=0.5)
            peakIFWHM2, indexes2 = beam.FWHM(t_grid, peakIPDF, frac=2)
            stdpeakIPDF = np.std(peakIPDF[indexes2])#(max(peakIPDF[indexes2]) - min(peakIPDF[indexes2]))/np.mean(peakIPDF[indexes2]) # Flat-top in the distribution!
            # print('stdpeakIPDF = ', stdpeakIPDF)
            # print 'Peak Fraction = ', 100*peakICDF[indexes][-1]-peakICDF[indexes][0], stdpeakIPDF
            beam.bin_time()
            t = 1e12*(beam.t - np.mean(beam.t))
            dt = beam.slice_length*(max(t) - min(t))

            sigmat = np.std(t)
            sigmap = np.std(beam.p)
            meanp = np.mean(beam.p)
            fitp = 100*sigmap/meanp
            peakI, peakIstd, peakIMomentumSpread, peakIEmittanceX, peakIEmittanceY, peakIMomentum, peakIDensity = beam.sliceAnalysis()
            peakI = max(peakIPDF)
            # chirp = beam.chirp
            chirp = 1e-6*(max(beam.cp) - min(beam.cp))

            print('bunchlength = ', 1e6*beam.sigma_z, 'um  chirp = ', -1*beam.chirp, 'MeV/ps')
            plt.subplots_adjust(left=0.1, bottom=0.1, right=0.9, top=0.9, wspace=0.25, hspace=0.15)

            ax1 = mw.getFigure().add_subplot(321)
            ax2 = mw.getFigure().add_subplot(323)
            ax3 = mw.getFigure().add_subplot(325)
            ax4 = mw.getFigure().add_subplot(322)
            ax5 = mw.getFigure().add_subplot(324)
            ax6 = mw.getFigure().add_subplot(326)

            ax = [[ax1,ax2,ax3],[ax4,ax5,ax6]]

            exponent = np.floor(np.log10(np.abs(beam.slice_length)))
            x = 10**(12) * np.array((beam.slice_bins - np.mean(beam.t)))
            ax[0][0].plot(x, beam.slice_peak_current)
            ax[0][0].set(xlabel='t (ps)', ylabel='I [A]')

            # ax[0][1].plot(t_grid, peakIPDF, color='blue', alpha=0.5, lw=3)
            # ax[0][1].fill_between(t_grid[indexes], peakIPDF[indexes], 0, facecolor='gray', edgecolor='gray', alpha=0.4)

            beam.read_HDF5_beam_file(dir+'/CLA-S07-CAV-04-SCR.hdf5')
            t = 1e12*(beam.t - np.mean(beam.t))
            dt = beam.slice_length*(max(t) - min(t))
            p = 1e-6*beam.cp
            ymax = max(p)+1
            ymin = min(p)-1
            if ymax - ymin < 5:
                ymax = np.mean(p) + 2.5
                ymin = np.mean(p) - 2.5
            ax[0][1].set_xlim(min(t) - dt, max(t) + dt)
            ax[0][1].hist2d(t, p, bins=(250,250), cmap=plt.cm.jet, range=[[min(t), max(t)],[ymin, ymax]])
            # ax[0][2].set_ylim(top=ymax, bottom=ymin)
            ax[0][1].set(ylabel='cp [Mev]')

            beam.read_HDF5_beam_file(dir+'/CLA-S07-APER-01.hdf5')
            t = 1e12*(beam.t - np.mean(beam.t))
            dt = beam.slice_length*(max(t) - min(t))
            p = 1e-6*beam.cp
            ymax = max(p)+1
            ymin = min(p)-1
            if ymax - ymin < 5:
                ymax = np.mean(p) + 2.5
                ymin = np.mean(p) - 2.5
            ax[0][2].set_xlim(min(t) - dt, max(t) + dt)
            ax[0][2].hist2d(t, p, bins=(250,250), cmap=plt.cm.jet, range=[[min(t), max(t)],[ymin, ymax]])
            # ax[0][2].set_ylim(top=ymax, bottom=ymin)
            ax[0][2].set(ylabel='cp [Mev]')

            ax[1][0].plot(x, 1e6*beam.slice_normalized_horizontal_emittance)
            ax[1][0].plot(x, 1e6*beam.slice_normalized_vertical_emittance)
            ax[1][0].set_ylim(top=3, bottom=0)
            ax[1][0].set(ylabel='emit_x / emit_y [m]')
            ax[1][0].grid()

            ax[1][1].plot(twiss.elegant['s'], 0.511*twiss.elegant['pCentral0'])
            # ax[1][1].set_ylim(top=1100, bottom=0)
            ax[1][1].set(ylabel='Momentum [MeV/c]')
            ax[1][1].grid()

            ax[1][2].plot(twiss.elegant['s'], 1e3*twiss['sigma_x'])
            ax[1][2].plot(twiss.elegant['s'], 1e3*twiss['sigma_y'])
            ax[1][2].set(ylabel='sigma_x / sigma_y [m]')
            ax[1][2].grid()

class OutLog:
    def __init__(self, edit, out=None, color=None):
        """(edit, out=None, color=None) -> can write stdout, stderr to a
        QTextEdit.
        edit = QTextEdit
        out = alternate stream ( can be the original sys.stdout )
        color = alternate color (i.e. color stderr a different color)
        """
        self.edit = edit
        self.out = None
        self.color = color

    def write(self, m):
        if self.color:
            tc = self.edit.textColor()
            self.edit.setTextColor(self.color)

        self.edit.moveCursor(qt.QTextCursor.End)
        self.edit.insertPlainText( m )

        if self.color:
            self.edit.setTextColor(tc)

        if self.out:
            self.out.write(m)

class mainApp(qt.QMainWindow):
    def __init__(self, parent = None):
        super(mainApp, self).__init__(parent)
        self.setMinimumWidth(1040/2)
        self.setMinimumHeight(1392/2)
        self.layout = qt.QGridLayout()
        self.widget = qt.QWidget()
        self.widget.setLayout(self.layout)
        self.setCentralWidget(self.widget)
        args = parser.parse_args()
        self.currentDirectory = os.path.abspath(args.directory)
        self.selectedDir = None

        self.table = qt.QListWidget()
        self.table.setMinimumWidth(200)
        self.table.setMaximumWidth(250)
        self.textedit = qt.QTextEdit()
        self.textedit.setMaximumHeight(150)
        self.textedit.setMinimumWidth(200)
        self.textedit.setMaximumWidth(250)
        self.update_directory_list()
        self.mw = MatplotlibWidget()
        self.mw.setMinimumWidth(1400)
        self.layout.addWidget(self.table,0,0,1,1)
        self.layout.addWidget(self.textedit,1,0,1,1)
        self.layout.addWidget(self.mw,0,1,2,3)

        # self.timer = qt.QTimer()
        # self.timer.timeout.connect(self.update_directory_list)
        # self.timer.start(10000)
        self.table.itemClicked.connect(self.analyse_image)
        self.table.itemDoubleClicked.connect(self.change_directory)

        # sys.stdout = OutLog( self.textedit, sys.stdout)
        # sys.stderr = OutLog( self.textedit, sys.stderr, qt.QColor(255,0,0) )

    def analyse_image(self, item):
        dir = self.currentDirectory + '/'+str(item.text())
        self.selectedDir = str(item.text())
        if not dir == '.' and not dir == '..':
            self.mw.getFigure().clear()
            # ''' initialise an instance of the stripPlot Widget '''
            analyse_image(self.mw, dir)
            self.mw.draw()

    def update_directory_list(self):
        index = None
        self.table.clear()
        start = time.time()
        # print (os.listdir(self.currentDirectory))
        dirs = [os.path.abspath(self.currentDirectory + '/' + a) for a in os.listdir(self.currentDirectory) if os.path.isdir(self.currentDirectory + '/' + a)]
        print ('read dirs in ', time.time() - start)
        if os.path.isfile(self.currentDirectory + '/best_solutions_running.csv'):
            solutions = []
            with open(self.currentDirectory + '/best_solutions_running.csv', 'rt') as csvfile:
                spamreader = csv.reader(csvfile, csv.QUOTE_NONE, delimiter=',')
                for row in spamreader:
                    solutions.append([float(a) for a in row])
            # solutions = np.array(sorted([s for s in solutions if s[-3] < 12 and s[-4] < 2], key=lambda a: a[-1]))
            solutions = np.array(sorted(solutions, key=lambda a: a[-1]))
            # print (solutions)
            iterdirs = [self.currentDirectory + '/iteration_' +str(int(a[-2])) for a in solutions if os.path.isdir(self.currentDirectory + '/iteration_' +str(int(a[-2]))) ]
            basedirs = [a for a in dirs if 'basefiles_' in a]
            vbcdirs = [a for a in dirs if 'vbc_' in a]
            setdirs = [a for a in dirs if 'set' in a]
            dirs = iterdirs# + setdirs + vbcdirs + basedirs
            # print 'sorted dirs in ', time.time() - start
        else:
            dirs.sort(key=os.path.getmtime, reverse=True)
        self.table.addItem('.')
        self.table.addItem('..')
        for i, d in enumerate(dirs[:100]):
            d = d.replace(self.currentDirectory+'/', '').replace(self.currentDirectory+'\\', '')
            item = self.table.addItem(d)
            if d == self.selectedDir:
                # print 'found dir = ', d, self.selectedDir
                index = i + 2
        if self.selectedDir is not None and index is not None:
            # try:
                self.table.itemClicked.disconnect(self.analyse_image)
                self.table.setCurrentRow(index)
                self.table.itemClicked.connect(self.analyse_image)
            # except:
            #     self.table.setCurrentRow(index)

    def change_directory(self, item):
        dir = str(item.text())
        print('changing directory! = ', os.path.abspath(self.currentDirectory + '/' + dir + '/'))
        self.currentDirectory = os.path.abspath(self.currentDirectory + '/' + dir + '/')
        self.update_directory_list()

def main():
    args = parser.parse_args()
    app = qt.QApplication(sys.argv)
    ex = mainApp()
    ex.show()
    sys.exit(app.exec_())

if __name__ == '__main__':
   main()
