import sys, os
import numpy as np
sys.path.append('./../../../')
from SimulationFramework.Examples.CLARA.Elegant.Optimise_longitudinal_Elegant import Optimise_Elegant
from SimulationFramework.Modules.merge_two_dicts import merge_two_dicts
from functools import partial
from ruamel import yaml
import SimulationFramework.Framework as fw

class FEBE(Optimise_Elegant):

    injector_startingvalues = [-9.,0.345,2.1e7,-16.,0.052500000000000005,-0.05]
    startingvalues = best = np.array([ 3.13845650e+07, -2.33062481e+01,  2.96752546e+07, -3.41502595e+00,
        2.74842883e+07,  1.87482967e+02,  3.11918859e+07,  5.07160187e+01,
       -1.22393267e-01,  6.12784140e-01])

    def __init__(self):
        super(FEBE, self).__init__()
        self.scaling = 6
        CLARA_dir = os.path.relpath(__file__+'/../../CLARA/')
        self.optfunc = partial(self.OptimisingFunction, scaling=self.scaling, post_injector=self.POST_INJECTOR, verbose=False,
            basefiles='../../../CLARA/basefiles_'+str(self.scaling)+'/', CLARA_dir=CLARA_dir, sample_interval=2**(3*1))
        self.verbose = True

    def load_best(self, filename):
        with open(filename, 'r') as infile:
            data = dict(yaml.load(infile, Loader=yaml.UnsafeLoader))
            # best = [data[n]['k1l'] for n in parameter_names]
            best = []
            for n, p in self.parameter_names:
                if n in data:
                    best.append(data[n][p])
                elif n == 'bunch_compressor' and p == 'set_angle':
                    best.append(data['CLA-VBC-MAG-DIP-01']['angle'])
                else:
                    print n, p
                    if not hasattr(self, 'framework'):
                        self.framework = fw.Framework(None)
                        self.framework.loadSettings(self.lattice)
                    best.append(self.framework[n][p])
            self.best = best

    def calculate_constraints(self):
        twiss = self.fit.twiss
        self.beam = self.fit.beam
        constraintsList = {}
        quadkls = self.fit.framework.getElementType('quadrupole','k1l')
        quadlengths = self.fit.framework.getElementType('quadrupole','length')
        quadnames = self.fit.framework.getElementType('quadrupole','objectname')

        self.beam.read_HDF5_beam_file(self.fit.dirname+'/CLA-FEB-W-FOCUS-01.hdf5')
        self.beam.slice_length = 0.025e-12
        self.beam.bin_time()
        sigmat = 1e12*np.std(self.beam.t)
        sigmap = np.std(self.beam.p)
        meanp = np.mean(self.beam.p)
        # emitx = 1e6*self.beam.normalized_horizontal_emittance
        # emity = 1e6*self.beam.normalized_horizontal_emittance
        # density = self.beam.density
        fitp = 100*sigmap/meanp
        peakI, peakIstd, peakIMomentumSpread, peakIEmittanceX, peakIEmittanceY, peakIMomentum, peakIDensity = self.beam.sliceAnalysis()

        linac_fields = np.array([1e-6*self.linac_fields[i] for i in [0,1,3]])

        twiss.read_elegant_twiss_files( self.fit.dirname+'/FEBE.twi' )
        ipindex = list(twiss.elegant['ElementName']).index('CLA-FEB-W-FOCUS-01')
        constraintsListFEBE = {
            # 'ip_enx': {'type': 'lessthan', 'value': 1e6*twiss.elegant['enx'][ipindex], 'limit': 2, 'weight': 0},
            # 'ip_eny': {'type': 'lessthan', 'value': 1e6*twiss.elegant['eny'][ipindex], 'limit': 0.5, 'weight': 2.5},
            'field_max': {'type': 'lessthan', 'value': linac_fields, 'limit': 32, 'weight': 30},
            'momentum': {'type': 'equalto', 'value': 0.511*twiss.elegant['pCentral0'][ipindex], 'limit': 250, 'weight': 1.5},
            'peakI_min': {'type': 'greaterthan', 'value': abs(peakI), 'limit': 1950, 'weight': 10},
            'peakI_max': {'type': 'lessthan', 'value': abs(peakI), 'limit': 2050, 'weight': 10},
            # 'peakIMomentumSpread': {'type': 'lessthan', 'value': peakIMomentumSpread, 'limit': 0.1, 'weight': 2},
            'peakIEmittanceX': {'type': 'lessthan', 'value': 1e6*peakIEmittanceX, 'limit': 2, 'weight': 1.5},
            'peakIEmittanceY': {'type': 'lessthan', 'value': 1e6*peakIEmittanceY, 'limit': 0.75, 'weight': 1.5},
            # 'peakIMomentum': {'type': 'equalto','value': 1e-6*peakIMomentum, 'limit': 250, 'weight': 20},
        }
        constraintsList = merge_two_dicts(constraintsList, constraintsListFEBE)

        # fitness = self.cons.constraints(constraintsList)
        if self.verbose:
            print self.cons.constraintsList(constraintsList)
        return constraintsList

opt = FEBE()
opt.set_changes_file(['./simplex_best_changes.yaml', './transverse_best_changes_upto_S07.yaml', './S07_transverse_best_changes.yaml', './FEBE_transverse_best_changes.yaml'])
opt.set_lattice_file('./FEBE.def')
opt.load_best('best_changes_1kA.yaml')
# opt.Nelder_Mead(step=[5e6, 10, 5e6, 10, 5e6, 10, 5e6, 10, 0.025, 0.5])
opt.Simplex()
