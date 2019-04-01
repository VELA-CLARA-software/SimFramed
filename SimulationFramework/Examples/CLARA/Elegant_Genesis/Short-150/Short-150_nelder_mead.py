import sys
import numpy as np
sys.path.append('./../')
from Optimise_Genesis_Elegant import Optimise_Genesis_Elegant

injector_startingvalues = [-9.,0.345,2.1e7,-16.,0.052500000000000005,-0.05]
startingvalues = best = np.array([17335579.78, -21.9040218, 15603832.7, -6.18551054, 13310042.9, 184.393421, 17989069.4, 44.386823, -0.128275132, 0.729080893])

def calculate_constraints(e, b, ee, be, l, g):
    brightness = (1e-4*e)/(1e-2*b)
    e = 1e2*e
    ee = 1e2*ee
    pos12m = list(g.z).index(10.)
    pos20m = list(g.z).index(20.)
    momentum = 1e-6*np.mean(g.momentum)
    print 'Momentum = ', momentum
    if e < 1:
        l = 500
    constraintsList = {
        'brightness': {'type': 'greaterthan', 'value': brightness, 'limit': 0.15, 'weight': 0},
        'bandwidth': {'type': 'lessthan', 'value': b, 'limit': 0.15, 'weight': 3},
        'pulse_energy': {'type': 'greaterthan', 'value': e, 'limit': 120, 'weight': 4},
        'bandwidth_end': {'type': 'lessthan', 'value': 1e2*abs(g.spectrum_lamdwidth_std[pos20m]), 'limit': 0.5, 'weight': 1},
        'pulse_energy_end': {'type': 'greaterthan', 'value': 1e6*abs(g.energy[pos20m]), 'limit': 200, 'weight': 2},
        'max_brightness_position': {'type': 'lessthan', 'value': abs(l), 'limit': 12, 'weight': 2.5},
        'min_brightness_position': {'type': 'greaterthan', 'value': abs(l), 'limit': 9, 'weight': 0.5},
        'energy_min': {'type': 'greaterthan', 'value': abs(momentum), 'limit': 145, 'weight': 50},
        'energy_max': {'type': 'lessthan', 'value': abs(momentum), 'limit': 155, 'weight': 50},
    }
    return constraintsList

opt = Optimise_Genesis_Elegant()
opt.calculate_constraints = calculate_constraints
opt.Nelder_Mead(best)
