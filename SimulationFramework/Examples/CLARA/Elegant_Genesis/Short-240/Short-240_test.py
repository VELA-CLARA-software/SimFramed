import sys, os
import numpy as np
sys.path.append('./../')
from Optimise_Genesis_Elegant import Optimise_Genesis_Elegant
from functools import partial

class Short240(Optimise_Genesis_Elegant):

    injector_startingvalues = [-9.,0.345,2.1e7,-16.,0.052500000000000005,-0.05]
    startingvalues = best = np.array([31635580, -21.9, 28403833, -6.2, 24210043, 184.4, 32589069, 44.4, -0.124, 0.93])
    # startingvalues = best = np.array([31295498.77911233, -21.796417988221847, 28533996.939464197, -4.3535424459076175, 25118020.516425118,
    # 186.03322780846452, 31693858.476702053, 49.0649623919314, -0.12195857313414843, 0.9046648051503599])
    startingvalues = best = np.array([ 3.13845650e+07, -2.33062481e+01,  2.96752546e+07, -3.41502595e+00,
        2.74842883e+07,  1.87482967e+02,  3.11918859e+07,  5.07160187e+01,
       -1.22393267e-01,  6.12784140e-01])

    def __init__(self):
        super(Short240, self).__init__()
        CLARA_dir = os.path.relpath(__file__+'/../../../../')
        scaling = 6
        self.optfunc = partial(self.OptimisingFunction, scaling=scaling, post_injector=self.POST_INJECTOR, verbose=False, basefiles='../../../../basefiles_'+str(scaling)+'/', CLARA_dir=CLARA_dir)

    def calculate_constraints(self):
        r = self.resultsDict
        pos12m = list(r['g'].z).index(10.)
        pos20m = list(r['g'].z).index(20.)
        print 'Momentum = ', r['momentum']
        constraintsList = {
            'brightness': {'type': 'greaterthan', 'value': r['brightness'], 'limit': 0.15, 'weight': 0},
            'bandwidth': {'type': 'lessthan', 'value': r['b'], 'limit': 0.075, 'weight': 3},
            'pulse_energy': {'type': 'greaterthan', 'value': 1e2*r['e'], 'limit': 150, 'weight': 4},
            'bandwidth_end': {'type': 'lessthan', 'value': 1e2*abs(r['g'].spectrum_lamdwidth_std[pos20m]), 'limit': 0.35, 'weight': 1},
            'pulse_energy_end': {'type': 'greaterthan', 'value': 1e6*abs(r['g'].energy[pos20m]), 'limit': 250, 'weight': 2},
            'max_brightness_position': {'type': 'lessthan', 'value': abs(r['l']), 'limit': 11, 'weight': 2.5},
            'min_brightness_position': {'type': 'greaterthan', 'value': abs(r['l']), 'limit': 9, 'weight': 0.5},
            'energy_min': {'type': 'greaterthan', 'value': abs(r['momentum']), 'limit': 240, 'weight': 50},
            'energy_max': {'type': 'lessthan', 'value': abs(r['momentum']), 'limit': 250, 'weight': 50},
            'linac_max': {'type': 'lessthan', 'value': 1e-6*self.linac_fields, 'limit': 32, 'weight': 2000},
        }
        return constraintsList

opt = Short240()
opt.parameter_names = [['CLA-S07-DCP-01', 'factor']]
opt.set_changes_file(['./nelder_mead_best_changes.yaml', './transverse_best_changes.yaml'])
opt.set_lattice_file('./Short-240.def')
opt.optdir = 'test/DCP_On'
opt.opt_iteration = ''
opt.optfunc([0.9]) 
#
# opt = Short240()
# opt.parameter_names = [['CLA-S07-DCP-01', 'factor']]
# opt.set_changes_file(['./nelder_mead_best_changes.yaml', './transverse_best_changes.yaml'])
# opt.set_lattice_file('./Short-240.def')
# opt.optdir = 'test/DCP_Off'
# opt.opt_iteration = ''
# opt.optfunc([0])
