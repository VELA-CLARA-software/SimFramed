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

def analysis(g):
    spectrum_lamdwidth_fwhm = np.zeros_like(g.z)
    spectrum_lamdwidth_std = np.zeros_like(g.z)
    for zz in range(g.nZ):
        spectrum_lamdwidth_fwhm[zz] = None
        spectrum_lamdwidth_std[zz] = None
        if np.sum(g.spec[:,zz])!=0:
            pos, width, arr = fwhm3(g.spec[:, zz])
            if width != None:
                if arr[0] == arr[-1]:
                    dlambda = abs(g.freq_lamd[pos] - g.freq_lamd[pos-1])
                else:
                    dlambda = abs( (g.freq_lamd[arr[0]] - g.freq_lamd[arr[-1]]) / (arr[0] - arr[-1]) )
                spectrum_lamdwidth_fwhm[zz] = dlambda * width / g.freq_lamd[pos]
                # spectrum_lamdwidth_fwhm[zz] = abs(g.freq_lamd[arr[0]] - g.freq_lamd[arr[-1]]) / g.freq_lamd[pos]  # the FWHM of spectral line (error when peakpos is at the edge of lamdscale)
            spectrum_lamdwidth_std[zz] = std_moment(g.freq_lamd, g.spec[:, zz]) / n_moment(g.freq_lamd, g.spec[:, zz], 0, 1)
    ######### end of section copied from genesis_plot.py
    brightness = g.energy / spectrum_lamdwidth_std
    for i in range(len(g.z)):
        if g.z[i] < 5:
            brightness[i] = 0
    return g.energy, spectrum_lamdwidth_std, g.z, brightness, g.xrms, g.yrms, g.I

def beam_analysis(d):
    beam = rbf.beam()
    beam.read_HDF5_beam_file(d+'/CLA-S07-APER-01.hdf5')
    return beam

def analyse_image(mw, dir=None):
    dir = str(dir)
    if dir is not None:
        if os.path.isfile(dir + '/' + 'run.0.gout'):
            args = parser.parse_args()
            out = read_out_file(dir + '/' + args.file, read_level=2)
            e, b, l, bright, rx, ry, current = analysis(out)
            mb = np.argmax(bright)
            print('bandwidth[br] = ', 1e2*b[mb], '%  pulse energy[br] =', 1e6*e[mb], 'uJ  Sat. Length[br] =', l[mb], 'm  Brightness[max] = ', bright[mb])
            en = np.argmax(e)
            print('bandwidth[en] = ', 1e2*b[en], '%  pulse energy[max] =', 1e6*e[en], 'uJ  Sat. Length[en] =', l[en], 'm  Brightness[en] = ', bright[en])
            beam = beam_analysis(dir)
            print('bunchlength = ', 1e6*beam.sigma_z, 'um  chirp = ', -1*beam.chirp, 'MeV/ps')
            plt.subplots_adjust(left=0.1, bottom=0.1, right=0.9, top=0.9, wspace=0.25, hspace=0.15)
            ax1 = mw.getFigure().add_subplot(421)
            ax1.plot(l, bright)
            ax1.set(ylabel='Brightness')
            ax1.grid()
            ax1.axvline(l[mb],color='black',ls='--')
            ax1.annotate(np.round(bright[mb], decimals=3), xy=(l[mb], bright[mb]), xytext=(3, 0.5*max(bright)), arrowprops=dict(facecolor='black', shrink=0.05))
            ax2 = mw.getFigure().add_subplot(423)
            ax2.plot(l, 1e6*e)
            ax2.set(ylabel='Pulse Energy [uJ]')
            ax2.grid()
            ax2.axvline(l[mb],color='black',ls='--')
            ax2.annotate(np.round(1e6*e[mb], decimals=3), xy=(l[mb], 1e6*e[mb]), xytext=(3, 0.5*max(1e6*e)), arrowprops=dict(facecolor='black', shrink=0.05))
            ax3 = mw.getFigure().add_subplot(425)
            ax3.plot(l, 1e2*b)
            ax3.set(xlabel='z (m)', ylabel='Bandwidth [%]')
            ax3.grid()
            ax3.axvline(l[mb],color='black',ls='--')
            ax3.annotate(np.round(1e2*b[mb], decimals=3), xy=(l[mb], 1e2*b[mb]), xytext=(3, 0.5*max(1e2*b)), arrowprops=dict(facecolor='black', shrink=0.05))
            # plt.show()
            ax6 = mw.getFigure().add_subplot(427)
            peakIpos = np.argmax(current)
            ax6.plot(l, 1e6*rx[peakIpos])
            ax6.plot(l, 1e6*ry[peakIpos])
            ax6.set(ylabel='x/y [um]')
            ax6.grid()

            ax4 = mw.getFigure().add_subplot(422)
            ax5 = mw.getFigure().add_subplot(424, sharex=ax4)
            ax4.get_shared_x_axes().join(ax4, ax5)
            t = 1e12*(beam.t - np.mean(beam.t))
            dt = 0.05*(max(t) - min(t))
            ax4.set_xlim(min(t) - dt, max(t) + dt)
            ax4.hist2d(1e12*(beam.t - np.mean(beam.t)), 1e-6*beam.cp, bins=(50,50), cmap=plt.cm.jet)
            # ax4.set_ylim(220,260)
            ax4.set(ylabel='cp [Mev]')

            # ax5.set_xlim(min(t) - dt, max(t) + dt)
            beam.slices = 100
            beam.bin_time()
            # beam.beam['total_charge'] = 2.5e-10
            print('charge = ', 1e12*beam.beam['total_charge'])
            exponent = np.floor(np.log10(np.abs(beam.slice_length)))
            x = 10**(12) * np.array((beam.slice_bins - np.mean(beam.t)))
            ax5.plot(x, beam.slice_peak_current)
            ax5.set(xlabel='t (ps)', ylabel='I [A]')

            ax7 = mw.getFigure().add_subplot(426, sharex=ax4)
            alphax=-0.189
            betax=3.76
            alphay=0
            betay=1.44
            # print 'beta_x = ', beam.beta_x, ' alpha_x = ', beam.alpha_x
            # print 'beta_y = ', beam.beta_y, ' alpha_y = ', beam.alpha_y
            beam.rematchXPlanePeakISlice(beta=betax, alpha=alphax)
            beam.rematchYPlanePeakISlice(beta=betay, alpha=alphay)
            # print 'beta_x = ', beam.beta_x, ' alpha_x = ', beam.alpha_x
            # print 'beta_y = ', beam.beta_y, ' alpha_y = ', beam.alpha_y
            ax7.plot(x, beam.slice_beta_x)
            ax7.plot(x, beam.slice_beta_y)
            ax7.set_ylim(top=20, bottom=0)
            ax7.set(ylabel='beta_x / beta_y [m]')
            ax7.grid()

            ax8 = mw.getFigure().add_subplot(428, sharex=ax4)
            ax8.plot(x, 1e6*beam.slice_normalized_horizontal_emittance)
            ax8.plot(x, 1e6*beam.slice_normalized_vertical_emittance)
            ax8.set(ylabel='enx / eny [um]')
            ax8.grid()

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

        self.timer = qt.QTimer()
        self.timer.timeout.connect(self.update_directory_list)
        self.timer.start(10000)
        self.table.itemClicked.connect(self.analyse_image)
        self.table.itemDoubleClicked.connect(self.change_directory)

        sys.stdout = OutLog( self.textedit, sys.stdout)
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
        dirs = [os.path.abspath(self.currentDirectory + '/' + a) for a in os.listdir(self.currentDirectory) if os.path.isdir(self.currentDirectory + '/' + a)]
        # print 'read dirs in ', time.time() - start
        if os.path.isfile(self.currentDirectory + '/best_solutions_running.csv'):
            solutions = []
            with open(self.currentDirectory + '/best_solutions_running.csv', 'rt') as csvfile:
                spamreader = csv.reader(csvfile, csv.QUOTE_NONE, delimiter=',')
                for row in spamreader:
                    solutions.append([float(a) for a in row])
            solutions = np.array(sorted(solutions, key=lambda a: a[-1]))
            iterdirs = [self.currentDirectory + '/iteration_' +str(int(a)) for a in solutions[:,-2] if os.path.isdir(self.currentDirectory + '/iteration_' +str(int(a))) ]
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
            try:
                self.table.itemClicked.disconnect(self.analyse_image)
                self.table.setCurrentRow(index)
                self.table.itemClicked.connect(self.analyse_image)
            except:
                self.table.setCurrentRow(index)

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
