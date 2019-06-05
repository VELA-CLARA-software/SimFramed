import os, sys
sys.path.append('../../../')
import numpy as np
import SimulationFramework.Framework as fw
from SimulationFramework.Modules.nelder_mead import nelder_mead
from SimulationFramework.Examples.CLARA.Elegant.Optimise_transverse import Optimise_transverse
from SimulationFramework.Modules.merge_two_dicts import merge_two_dicts
from ruamel import yaml
import matplotlib.pyplot as plt
plt.ion()
plt.show()
fig, ax = plt.subplots(2, 3, sharey=False,
                       figsize=(16, 6))

framework = fw.Framework(None)
framework.loadSettings('Lattices/claraX400_v12_80MVm_Elegant.def')
parameters = framework.getElementType('quadrupole','k1l')
names = framework.getElementType('quadrupole','objectname')
index1 = 0# names.index('CLA-S06-MAG-QUAD-01')
parameter_names = []
parameter_names += [q for q in names[index1:]]
best = parameters

with open('transverse_best_changes.yaml', 'r') as infile:
    data = dict(yaml.load(infile, Loader=yaml.UnsafeLoader))
    # best = [data[n]['k1l'] for n in parameter_names]
    best = []
    for n in parameter_names:
        if n in data:
            best.append(data[n]['k1l'])
        else:
            print(n)
            best.append(framework[n]['k1l'])

class XARA_Transverse(Optimise_transverse):

    def __init__(self, lattice='Lattices/claraX400_v12_80MVm.def', scaling=6):
        super(XARA_Transverse, self).__init__(lattice=lattice, scaling=scaling)
        names = framework.getElementType('quadrupole','objectname')
        # self.framework.change_Lattice_Code('All','elegant')
        self.scaling = scaling
        self.parameter_names = []
        self.parameters = []
        self.parameter_names += [q for q in names[index1:]]
        self.parameters += [[q, 'k1l'] for q in names[index1:]]
        self.save_parameters = []
        self.CLARA_dir = os.path.relpath(__file__+'/../../../../CLARA/')
        # self.base_files = '../../CLARA/basefiles_' + str(int(scaling)) + '/'
        self.best_changes = './transverse_best_changes.yaml'
        self.start_file = 'CLARAX'

    def calculateBeamParameters(self):
        twiss = self.twiss
        self.framework[self.start_file].prefix = self.CLARA_dir+'/basefiles_'+str(self.scaling)+'/'
        # print('prefix basedir = ', self.CLARA_dir, self.framework[self.start_file].prefix)
        self.framework[self.start_file].sample_interval = 2**(3*2)
        # print ('###############  STARTING TRACKING  ###############')
        self.framework.track(startfile=self.start_file)
        # print ('###############  ENDING TRACKING  ###############')

        constraintsList = {}

        twiss.reset_dicts()

        for lat in ['CLARAX']:
            quadkls = self.framework[lat].getElementType('quadrupole','k1l')
            quadlengths = self.framework[lat].getElementType('quadrupole','length')
            constraintsListQuads = {
                'max_k_'+lat: {'type': 'lessthan', 'value': [abs(k) for k, l in zip(quadkls, quadlengths)], 'limit': 2.0, 'weight': 75},

            }
            constraintsList = merge_two_dicts(constraintsList, constraintsListQuads)

        twiss.read_elegant_twiss_files( [ self.dirname+'/CLARAX.twi' ])
        ipindex1 = list(twiss['element_name']).index('CLA-S05-MARK-01')
        ipindex2 = list(twiss['element_name']).index('CLA-S06-MARK-01')
        constraintsListFEBE = {
            'min_sigmax_1': {'type': 'greaterthan', 'value': 1e3*twiss['sigma_x'][:ipindex1], 'limit': 0.25, 'weight': 15},
            'min_sigmax_2': {'type': 'greaterthan', 'value': 1e3*twiss['sigma_x'][ipindex2:], 'limit': 0.05, 'weight': 15},
            'min_sigmay_1': {'type': 'greaterthan', 'value': 1e3*twiss['sigma_y'][:ipindex1], 'limit': 0.25, 'weight': 15},
            'min_sigmay_2': {'type': 'greaterthan', 'value': 1e3*twiss['sigma_y'][ipindex2:], 'limit': 0.05, 'weight': 15},

            'max_sigmax_1': {'type': 'lessthan', 'value': 1e3*twiss['sigma_x'][:ipindex1], 'limit': 1, 'weight': 15},
            'max_sigmax_2': {'type': 'lessthan', 'value': 1e3*twiss['sigma_x'][ipindex2:], 'limit': 0.25, 'weight': 15},
            'max_sigmay_1': {'type': 'lessthan', 'value': 1e3*twiss['sigma_y'][:ipindex1], 'limit': 1, 'weight': 15},
            'max_sigmay_2': {'type': 'lessthan', 'value': 1e3*twiss['sigma_y'][ipindex2:], 'limit': 0.25, 'weight': 15},

            'max_emitx_2': {'type': 'lessthan', 'value': 1e6*twiss['enx'][ipindex2:], 'limit': 2.0, 'weight': 0},
            'max_emity_2': {'type': 'lessthan', 'value': 1e6*twiss['eny'][ipindex2:], 'limit': 0.6, 'weight': 15},
        }
        constraintsList = merge_two_dicts(constraintsList, constraintsListFEBE)

        self.beam.read_HDF5_beam_file(self.dirname+'/CLA-S07-APER-01.hdf5')
        self.beam.slice_length = 0.01e-12

        t = 1e12*(self.beam.t-np.mean(self.beam.t))
        t_grid = np.linspace(min(t), max(t), 2**8)
        peakIPDF = self.beam.PDF(t, t_grid, bandwidth=self.beam.rms(t)/(2**4))
        peakICDF = self.beam.CDF(t, t_grid, bandwidth=self.beam.rms(t)/(2**4))
        peakIFWHM, indexes = self.beam.FWHM(t_grid, peakIPDF, frac=0.5)
        peakIFWHM2, indexes2 = self.beam.FWHM(t_grid, peakIPDF, frac=2)
        stdpeakIPDF = (max(peakIPDF[indexes2]) - min(peakIPDF[indexes2]))/np.mean(peakIPDF[indexes2])
        # print('stdpeakIPDF = ', stdpeakIPDF)
        # print 'Peak Fraction = ', 100*peakICDF[indexes][-1]-peakICDF[indexes][0], stdpeakIPDF
        self.beam.bin_time()
        t = 1e12*(self.beam.t - np.mean(self.beam.t))
        dt = 0.05*(max(t) - min(t))

        self.constraintsList = constraintsList
        fitness = self.cons.constraints(constraintsList)
        if self.verbose:
            print(self.cons.constraintsList(constraintsList))

        if fitness < self.bestfit:

            ax[0][0].clear()
            ax[0][1].clear()
            ax[0][2].clear()
            ax[1][0].clear()
            ax[1][1].clear()
            ax[1][2].clear()

            exponent = np.floor(np.log10(np.abs(self.beam.slice_length)))
            x = 10**(12) * np.array((self.beam.slice_bins - np.mean(self.beam.t)))
            ax[0][0].plot(x, self.beam.slice_peak_current)
            ax[0][0].set(xlabel='t (ps)', ylabel='I [A]')

            ax[0][1].plot(t_grid, peakIPDF, color='blue', alpha=0.5, lw=3)
            ax[0][1].fill_between(t_grid[indexes], peakIPDF[indexes], 0, facecolor='gray', edgecolor='gray', alpha=0.4)

            # ax[0][1].plot(t_grid, peakICDF, color='blue', alpha=0.5, lw=3)
            # ax[0][1].fill_between(t_grid[indexes], peakICDF[indexes], 0, facecolor='gray', edgecolor='gray', alpha=0.4)

            ax[0][2].set_xlim(min(t) - dt, max(t) + dt)
            t = 1e12*(self.beam.t - np.mean(self.beam.t))
            p = 1e-6*self.beam.cp
            ymax = max(p)+1
            ymin = min(p)-1
            if ymax - ymin < 5:
                ymax = np.mean(p) + 2.5
                ymin = np.mean(p) - 2.5
            ax[0][2].hist2d(t, p, bins=(50,50), cmap=plt.cm.jet, range=[[min(t), max(t)],[ymin, ymax]])
            # ax[0][2].set_ylim(top=ymax, bottom=ymin)
            ax[0][2].set(ylabel='cp [Mev]')

            ax[1][0].plot(twiss.elegant['s'], 1e6*twiss.elegant['enx'])
            ax[1][0].plot(twiss.elegant['s'], 1e6*twiss.elegant['eny'])
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

        plt.draw()
        plt.pause(0.001)

        return fitness


if __name__ == "__main__":
        fit = XARA_Transverse(lattice='Lattices/claraX400_v12_80MVm_Elegant.def', scaling=6)
        fit.setChangesFile(['./nelder_mead_best_changes.yaml'])
        fit.verbose = False
        fit.Nelder_Mead(best, step=0.01)
        # fit.Simplex(best)
