import sys, os
import numpy as np
sys.path.append('./../../../')
from SimulationFramework.ClassFiles.Optimise_longitudinal_Elegant import Optimise_Elegant
from SimulationFramework.Modules.merge_two_dicts import merge_two_dicts
from functools import partial
from ruamel import yaml
import SimulationFramework.Framework as fw
import matplotlib.pyplot as plt
import argparse

parser = argparse.ArgumentParser(description='Perform an optimisation.')
parser.add_argument('type', default='nelder_mead')

class FEBE(Optimise_Elegant):

    # injector_startingvalues = [-9.,0.345,2.1e7,-16.,0.052500000000000005,-0.05]
    # startingvalues = best = np.array([ 3.13845650e+07, -2.33062481e+01,  2.96752546e+07, -3.41502595e+00,
    #     2.74842883e+07,  1.87482967e+02,  3.11918859e+07,  5.07160187e+01,
    #    -1.22393267e-01,  6.12784140e-01])

    parameter_names = [
        # ['startcharge','charge'],
        ['CLA-L01-CAV', 'field_amplitude'],
        ['CLA-L01-CAV', 'phase'],
        ['CLA-L02-CAV', 'field_amplitude'],
        ['CLA-L02-CAV', 'phase'],
        ['CLA-L03-CAV', 'field_amplitude'],
        ['CLA-L03-CAV', 'phase'],
        ['CLA-L4H-CAV', 'field_amplitude'],
        ['CLA-L4H-CAV', 'phase'],
        ['CLA-L04-CAV', 'field_amplitude'],
        ['CLA-L04-CAV', 'phase'],
        ['bunch_compressor', 'angle'],
        ['CLA-S07-DCP-01', 'factor'],
        # ['FODO_D', 'k1l'],
        # ['FODO_F', 'k1l'],
    ]

    def __init__(self):
        super(FEBE, self).__init__()
        self.parameter_names.append(['FODO_F', 'k1l'])
        self.parameter_names.append(['FODO_D', 'k1l'])
        self.scaling = 6
        self.sample_interval = 2**(3*0)
        self.base_files = '../../basefiles_'+str(self.scaling)+'/'
        self.clean = True
        self.doTracking = True
        self.startcharge = 250

    def before_tracking(self):
            self.framework.defineElegantCommand(scaling=6)
            csrbins = int(round(2**(3*self.scaling) / self.sample_interval / 128, -1))
            lscbins = int(round(2**(3*self.scaling) / self.sample_interval / 256, -1))
            elements = self.framework.elementObjects.values()
            for e in elements:
                e.lsc_enable = True
                e.lsc_bins = lscbins
                e.current_bins = 0
                e.csr_bins = csrbins
                e.longitudinal_wakefield_enable = True
                e.transverse_wakefield_enable = True
                e.smoothing_half_width = 1
                e.lsc_high_frequency_cutoff_start = 0.2
                e.lsc_high_frequency_cutoff_end = 0.25
                e.smoothing = 1
                # e.end1_focus = 0
                # e.end2_focus = 0
                # e.body_focus_model = "None"
            lattices = self.framework.latticeObjects.values()
            for l in lattices:
                l.lscDrifts = True
                l.lsc_bins = lscbins
                l.lsc_high_frequency_cutoff_start = 0.2
                l.lsc_high_frequency_cutoff_end = 0.25
                l.smoothing_half_width = 1
                l.smoothing = 1
            # self.framework['FEBE'].betax = 0.748377
            # self.framework['FEBE'].betay = 3.96107
            # self.framework['FEBE'].alphax = -0.61447
            # self.framework['FEBE'].alphay = 0.872954

    def calculate_constraints(self):
        constraintsList = {}
        quadkls = self.framework.getElementType('quadrupole','k1l')
        quadlengths = self.framework.getElementType('quadrupole','length')
        quadnames = self.framework.getElementType('quadrupole','objectname')

        self.beam.read_HDF5_beam_file(self.dirname+'/CLA-FEH-FOCUS-01.hdf5')
        self.beam.slice_length = 0.05e-12

        t = 1e12*(self.beam.t-np.mean(self.beam.t))
        t_grid = np.linspace(min(t), max(t), 2**8)
        peakIPDF = self.beam.PDF(t, t_grid, bandwidth=self.beam.rms(t)/(2**4))
        peakICDF = self.beam.CDF(t, t_grid, bandwidth=self.beam.rms(t)/(2**4))
        peakIFWHM, indexes = self.beam.FWHM(t_grid, peakIPDF, frac=4)
        # peakIFWHM2, indexes2 = self.beam.FWHM(t_grid, peakIPDF, frac=10)
        # stdpeakIPDF = (max(peakIPDF[indexes2]) - min(peakIPDF[indexes2]))/np.mean(peakIPDF[indexes2])
        # print('stdpeakIPDF = ', stdpeakIPDF)
        # print 'Peak Fraction = ', 100*peakICDF[indexes][-1]-peakICDF[indexes][0], peakICDF[indexes][-1], peakICDF[indexes][0]
        #
        # import matplotlib.pyplot as plt
        # fig, ax = plt.subplots(1, 2, sharey=False,
        #                        figsize=(13, 3))
        #
        # i=0
        # ax[0].plot(t_grid, peakIPDF, color='blue', alpha=0.5, lw=3)
        # ax[0].fill_between(t_grid[indexes], peakIPDF[indexes], 0, facecolor='gray', edgecolor='gray', alpha=0.4)
        # ax[1].plot(t_grid, peakICDF, color='blue', alpha=0.5, lw=3)
        # ax[1].fill_between(t_grid[indexes], peakICDF[indexes], 0, facecolor='gray', edgecolor='gray', alpha=0.4)
        # plt.show()
        # exit()
        self.beam.bin_time()
        sigmat = np.std(t)
        sigmap = np.std(self.beam.p)
        meanp = np.mean(self.beam.p)
        # emitx = 1e6*self.beam.normalized_horizontal_emittance
        # emity = 1e6*self.beam.normalized_horizontal_emittance
        # density = self.beam.density
        fitp = 100*sigmap/meanp
        peakI, peakIstd, peakIMomentumSpread, peakIEmittanceX, peakIEmittanceY, peakIMomentum, peakIDensity = self.beam.sliceAnalysis()

        fhc_field = np.array([1e-6*self.linac_fields[i] for i in [3]])
        linac_fields = np.array([1e-6*self.linac_fields[i] for i in [0,1,2, 4]])

        self.twiss.read_elegant_twiss_files( self.dirname+'/FEBE.twi' )
        ipindex = list(self.twiss.elegant['ElementName']).index('CLA-FEH-FOCUS-01')
        ipindex2 = list(self.twiss.elegant['ElementName']).index('CLA-FEH-FOCUS-02')
        constraintsListFEBE = {
            # 'ip_enx': {'type': 'lessthan', 'value': 1e6*self.twiss.elegant['enx'][ipindex], 'limit': 2, 'weight': 0},
            # 'ip_eny': {'type': 'lessthan', 'value': 1e6*self.twiss.elegant['eny'][ipindex], 'limit': 0.5, 'weight': 2.5},
            'field_max': {'type': 'lessthan', 'value': linac_fields, 'limit': 32, 'weight': 3000},
            'field_max_fhc': {'type': 'lessthan', 'value': fhc_field, 'limit': 50, 'weight': 3000},
            'momentum_max': {'type': 'lessthan', 'value': 0.511*self.twiss.elegant['pCentral0'][ipindex], 'limit': 260, 'weight': 250},
            'momentum_min': {'type': 'greaterthan', 'value': 0.511*self.twiss.elegant['pCentral0'][ipindex], 'limit': 240, 'weight': 150},
            # 'sigma_x_IP': {'type': 'lessthan', 'value': self.twiss.elegant['Sx'][ipindex], 'limit': 15e-6, 'weight': 15},
            # 'sigma_y_IP': {'type': 'lessthan', 'value': self.twiss.elegant['Sy'][ipindex], 'limit': 8e-6, 'weight': 15},
            'sigma_x_max': {'type': 'lessthan', 'value': 1e3*6*self.twiss.elegant['Sx'][:ipindex2], 'limit': 37, 'weight': 150},
            'sigma_y_max': {'type': 'lessthan', 'value': 1e3*6*self.twiss.elegant['Sy'][:ipindex2], 'limit': 37, 'weight': 150},
            'peakI_min': {'type': 'greaterthan', 'value': abs(peakI), 'limit': 5000, 'weight': 100},
            # 'peakI_max': {'type': 'lessthan', 'value': abs(peakI), 'limit': 1050, 'weight': 20},
            # 'peakIMomentumSpread': {'type': 'lessthan', 'value': peakIMomentumSpread, 'limit': 0.1, 'weight': 2},
            'peakIEmittanceX': {'type': 'lessthan', 'value': 1e6*peakIEmittanceX, 'limit': 25, 'weight': 1.5},
            'peakIEmittanceY': {'type': 'lessthan', 'value': 1e6*peakIEmittanceY, 'limit': 1, 'weight': 1.5},
            'peakIFWHM': {'type': 'lessthan','value': peakIFWHM, 'limit': 0.03, 'weight': 50},
            # 'stdpeakIFWHM': {'type': 'lessthan','value': stdpeakIPDF, 'limit': 1, 'weight': 0},
            'peakIFraction': {'type': 'greaterthan','value': 100*peakICDF[indexes][-1]-peakICDF[indexes][0], 'limit': 75, 'weight': 25},
        }
        constraintsList = merge_two_dicts(constraintsList, constraintsListFEBE)

        # fitness = self.cons.constraints(constraintsList)
        if self.verbose:
            print(self.cons.constraintsList(constraintsList))
        return constraintsList

if __name__ == "__main__":
    opt = FEBE()
    opt.set_changes_file(['./transverse_best_changes_upto_S07_250pC.yaml', './S07_transverse_best_changes_250pC.yaml', './FEBE_transverse_best_changes.yaml'])
    opt.set_lattice_file('./FEBE_Single_L01.def')
    opt.set_start_file('PreFEBE')
    args = parser.parse_args()
    if args.type == 'nelder_mead':
        print('Performing a Nelder-Mead Optimisation...')
        opt.load_best('./nelder_mead_best_changes_250pC.yaml')
        opt.Nelder_Mead(best_changes='./nelder_mead_best_changes_250pC.yaml', subdir='nelder_mead_250pC', step=[1e6, 2, 1e6, 2, 1e6, 2, 1e6, 2, 1e6, 2, 0.001, 0.02, 0.1,0.1])
    elif args.type == 'simplex':
        print('Performing a SciPy Simplex Optimisation...')
        opt.load_best('./simplex_best_changes_250pC.yaml')
        opt.Simplex(best_changes='./simplex_best_changes_250pC.yaml', subdir='simplex_250pC', maxiter=300)
