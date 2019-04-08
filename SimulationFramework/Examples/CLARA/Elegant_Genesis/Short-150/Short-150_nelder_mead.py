import sys
import numpy as np
sys.path.append('./../')
from Optimise_Genesis_Elegant import Optimise_Genesis_Elegant


class Short150(Optimise_Genesis_Elegant):

    injector_startingvalues = [-9.,0.345,2.1e7,-16.,0.052500000000000005,-0.05]
    startingvalues = best = np.array([17582732,-24.7,15749127,-2.8,12816219,184.7,17950455,49.3,-0.1232,0.375])
    startingvalues = best = [ 1.76313906e+07, -2.35583473e+01,  1.52605218e+07, -1.23994485e+00,
        1.30308938e+07,  1.85760004e+02,  1.85635741e+07,  5.35698440e+01,
       -1.22352171e-01,  3.45711353e-01]

    def calculate_constraints(self):
        r = self.resultsDict
        pos12m = list(r['g'].z).index(10.)
        pos20m = list(r['g'].z).index(20.)
        print 'Momentum = ', r['momentum']
        constraintsList = {
            'brightness': {'type': 'greaterthan', 'value': r['brightness'], 'limit': 0.1, 'weight': 0},
            'bandwidth': {'type': 'lessthan', 'value': r['b'], 'limit': 0.2, 'weight': 3},
            'pulse_energy': {'type': 'greaterthan', 'value': 1e2*r['e'], 'limit': 120, 'weight': 4},
            'bandwidth_end': {'type': 'lessthan', 'value': 1e2*abs(r['g'].spectrum_lamdwidth_std[pos20m]), 'limit': 0.5, 'weight': 1},
            'pulse_energy_end': {'type': 'greaterthan', 'value': 1e6*abs(r['g'].energy[pos20m]), 'limit': 150, 'weight': 2},
            'max_brightness_position': {'type': 'lessthan', 'value': abs(r['l']), 'limit': 8, 'weight': 2.5},
            'min_brightness_position': {'type': 'greaterthan', 'value': abs(r['l']), 'limit': 6, 'weight': 0.5},
            'energy_min': {'type': 'greaterthan', 'value': abs(r['momentum']), 'limit': 145, 'weight': 50},
            'energy_max': {'type': 'lessthan', 'value': abs(r['momentum']), 'limit': 155, 'weight': 50},
            'linac_max': {'type': 'lessthan', 'value': 1e-6*self.linac_fields, 'limit': 32, 'weight': 2000},
        }
        return constraintsList

opt = Short150()
opt.set_changes_file('transverse_best_changes.yaml')
opt.Nelder_Mead(step=[1e6, 2, 1e6, 2, 1e6, 2, 1e6, 2, 0.01, 0.1])
