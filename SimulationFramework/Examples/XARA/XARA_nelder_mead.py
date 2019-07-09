import sys, os
import numpy as np
sys.path.append('./../../../')
from SimulationFramework.ClassFiles.Optimise_longitudinal_Elegant import Optimise_Elegant, saveState
from SimulationFramework.Modules.merge_two_dicts import merge_two_dicts
from functools import partial
from ruamel import yaml
import SimulationFramework.Framework as fw

import matplotlib.pyplot as plt

class XARA(Optimise_Elegant):

    parameter_names = [
        # ['CLA-L01-CAV', 'field_amplitude'],
        # ['CLA-L01-CAV', 'phase'],
        ['CLA-L02-CAV', 'field_amplitude'],
        ['CLA-L02-CAV', 'phase'],
        ['CLA-L03-CAV', 'field_amplitude'],
        ['CLA-L03-CAV', 'phase'],
        ['CLA-L4H-CAV', 'field_amplitude'],
        ['CLA-L4H-CAV', 'phase'],
        ['CLA-L04-CAV-01', 'field_amplitude'],
        ['CLA-L04-CAV-01', 'phase'],
        ['CLA-L04-CAV-02', 'field_amplitude'],
        ['CLA-L04-CAV-02', 'phase'],
        ['CLA-L04-CAV-03', 'field_amplitude'],
        ['CLA-L04-CAV-03', 'phase'],
        ['CLA-L04-CAV-04', 'field_amplitude'],
        ['CLA-L04-CAV-04', 'phase'],
        ['bunch_compressor', 'angle'],
        ['CLA-S07-DCP-01', 'factor'],
        ['CLA-S04-LH-SCA', 'relative_momentum_scatter'],
    ]

    def __init__(self):
        super(XARA, self).__init__()
        self.scaling = 6
        self.CLARA_dir = os.path.relpath(__file__+'/../../../../CLARA/')
        self.sample_interval=2**(3*2)
        # self.optfunc = partial(self.OptimisingFunction, scaling=self.scaling, post_injector=self.POST_INJECTOR, verbose=False, sample_interval=2**(3*2))
        self.verbose = False
        self.plotting = False

    def calculate_constraints(self):
        beam = self.beam
        constraintsList = {}
        quadkls = self.framework.getElementType('quadrupole','k1l')
        quadlengths = self.framework.getElementType('quadrupole','length')
        quadnames = self.framework.getElementType('quadrupole','objectname')

        self.beam.read_HDF5_beam_file(self.dirname+'/CLA-S07-APER-01.hdf5')
        self.beam.slice_length = 0.01e-12

        t = 1e12*(self.beam.t-np.mean(self.beam.t))
        t_grid = np.linspace(min(t), max(t), 2**8)
        # t_grid = np.arange(min(t), max(t), 0.01)
        peakIPDF = self.beam.PDFI(t, t_grid, bandwidth=self.beam.rms(t)/(2**3))*250
        peakICDF = self.beam.CDF(t, t_grid, bandwidth=self.beam.rms(t)/(2**3))
        peakIFWHM, indexes = self.beam.FWHM(t_grid, peakIPDF, frac=0.5)
        peakIFWHM2, indexes2 = self.beam.FWHM(t_grid, peakIPDF, frac=2)
        stdpeakIPDF = np.std(peakIPDF[indexes2])#(max(peakIPDF[indexes2]) - min(peakIPDF[indexes2]))/np.mean(peakIPDF[indexes2]) # Flat-top in the distribution!
        # print('stdpeakIPDF = ', stdpeakIPDF)
        # print 'Peak Fraction = ', 100*peakICDF[indexes][-1]-peakICDF[indexes][0], stdpeakIPDF
        self.beam.bin_time()
        t = 1e12*(self.beam.t - np.mean(self.beam.t))
        dt = self.beam.slice_length*(max(t) - min(t))

        sigmat = np.std(t)
        sigmap = np.std(self.beam.p)
        meanp = np.mean(self.beam.p)
        fitp = 100*sigmap/meanp
        peakI, peakIstd, peakIMomentumSpread, peakIEmittanceX, peakIEmittanceY, peakIMomentum, peakIDensity = self.beam.sliceAnalysis()
        peakI = max(peakIPDF)
        # chirp = self.beam.chirp
        chirp = 1e-6*(max(self.beam.cp) - min(self.beam.cp))

        slinac_fields = np.array([1e-6*self.linac_fields[i] for i in [0,1]])
        x4hlinac_fields = np.array([1e-6*self.linac_fields[i] for i in [2]])
        xlinac_fields = np.array([1e-6*self.linac_fields[i] for i in [4,5,6]])

        self.twiss.read_elegant_twiss_files( self.dirname+'/CLARAX.twi' )
        ipindex = list(self.twiss.elegant['ElementName']).index('CLA-S07-APER-01')
        constraintsListXARA = {
            # 'ip_enx': {'type': 'lessthan', 'value': 1e6*self.twiss.elegant['enx'][ipindex], 'limit': 2, 'weight': 0},
            # 'ip_eny': {'type': 'lessthan', 'value': 1e6*self.twiss.elegant['eny'][ipindex], 'limit': 0.5, 'weight': 2.5},
            'field_max_s': {'type': 'lessthan', 'value': slinac_fields, 'limit': 32, 'weight': 500},
            'field_max_x': {'type': 'lessthan', 'value': xlinac_fields, 'limit': 90, 'weight': 500},
            'field_max_x4h': {'type': 'lessthan', 'value': x4hlinac_fields, 'limit': 35, 'weight': 500},
            'momentum_max': {'type': 'lessthan', 'value': 0.511*self.twiss.elegant['pCentral0'][ipindex], 'limit': 1050, 'weight': 250},
            'momentum_min': {'type': 'greaterthan', 'value': 0.511*self.twiss.elegant['pCentral0'][ipindex], 'limit': 990, 'weight': 150},
            'peakI_min': {'type': 'greaterthan', 'value': abs(peakI), 'limit': 1000, 'weight': 2000},
            'peakI_max': {'type': 'lessthan', 'value': abs(peakI), 'limit': 1025, 'weight': 2000},
            # 'peakIMomentumSpread': {'type': 'lessthan', 'value': peakIMomentumSpread, 'limit': 0.1, 'weight': 2},
            'peakIEmittanceX': {'type': 'lessthan', 'value': 1e6*peakIEmittanceX, 'limit': 0.75, 'weight': 15},
            'peakIEmittanceY': {'type': 'lessthan', 'value': 1e6*peakIEmittanceY, 'limit': 0.75, 'weight': 1.5},
            'peakIFWHM': {'type': 'lessthan','value': peakIFWHM, 'limit': 1, 'weight': 100},
            # 'stdpeakIFWHM': {'type': 'lessthan','value': stdpeakIPDF, 'limit': 30, 'weight': 50},
            'peakIFraction': {'type': 'greaterthan','value': 100*peakICDF[indexes][-1]-peakICDF[indexes][0], 'limit': 90, 'weight': 20},
            'chirp': {'type': 'lessthan', 'value': abs(chirp), 'limit': 0.5, 'weight': 25},
        }
        constraintsList = merge_two_dicts(constraintsList, constraintsListXARA)

        fitness = self.cons.constraints(constraintsList)

        if self.plotting and fitness < self.bestfit:

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
            ax[0][2].hist2d(t, p, bins=(250,250), cmap=plt.cm.jet, range=[[min(t), max(t)],[ymin, ymax]])
            # ax[0][2].set_ylim(top=ymax, bottom=ymin)
            ax[0][2].set(ylabel='cp [Mev]')

            ax[1][0].plot(x, 1e6*self.beam.slice_normalized_horizontal_emittance)
            ax[1][0].plot(x, 1e6*self.beam.slice_normalized_vertical_emittance)
            ax[1][0].set_ylim(top=3, bottom=0)
            ax[1][0].set(ylabel='emit_x / emit_y [m]')
            ax[1][0].grid()

            ax[1][1].plot(self.twiss.elegant['s'], 0.511*self.twiss.elegant['pCentral0'])
            # ax[1][1].set_ylim(top=1100, bottom=0)
            ax[1][1].set(ylabel='Momentum [MeV/c]')
            ax[1][1].grid()

            ax[1][2].plot(self.twiss.elegant['s'], 1e3*self.twiss['sigma_x'])
            ax[1][2].plot(self.twiss.elegant['s'], 1e3*self.twiss['sigma_y'])
            ax[1][2].set(ylabel='sigma_x / sigma_y [m]')
            ax[1][2].grid()

        # fig.canvas.draw_idle()
        # plt.draw()
        if self.plotting:
            plt.pause(0.001)

        if self.verbose:
            print(self.cons.constraintsList(constraintsList))
        return constraintsList

    def before_tracking(self):
        elements = self.framework.elementObjects.values()
        for e in elements:
            e.lsc_enable = True
            e.lsc_bins = 100
            e.current_bins = 0
            e.longitudinal_wakefield_enable = True
            e.transverse_wakefield_enable = True
        lattices = self.framework.latticeObjects.values()
        for l in lattices:
            l.lscDrifts = True
            l.lsc_bins = 100

if __name__ == "__main__":
    plt.ion()
    plt.show()
    fig, ax = plt.subplots(2, 3, sharey=False,
                           figsize=(16, 6))

    opt = XARA()
    opt.set_changes_file(['nelder_mead_best_changes.yaml', './transverse_best_changes.yaml'])
    opt.set_lattice_file('Lattices/claraX400_v12_80MVm_Elegant.def')
    opt.set_start_file('CLARAX')
    opt.load_best('nelder_mead_best_changes.yaml')
    opt.plotting = True
    opt.Nelder_Mead(step=[5e6, 5, 5e6, 5, 5e6, 5, 5e6, 5, 5e6, 5,  5e6, 5,  5e6, 5, 0.005, 0.1, 5e-5])
    # opt.Simplex()
