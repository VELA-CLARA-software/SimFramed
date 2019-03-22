import os, errno, sys
import numpy as np
import random
sys.path.append(os.path.abspath(__file__+'/../../../../../'))
from SimulationFramework.Modules.constraints import constraintsClass
import time
import csv
import copy
import genesisBeamFile
from functools import partial
from scipy.optimize import minimize

def saveState(args, n, fitness):
    with open('simplex/best_solutions_running.csv','a') as out:
        csv_out=csv.writer(out)
        args=list(args)
        args.append(n)
        args.append(fitness)
        csv_out.writerow(args)

def saveParameterFile(best, file='clara_elegantgenesis_best.yaml'):
    allparams = zip(*(injparameternames+parameternames))
    output = {}
    for p, k, v in zip(allparams[0], allparams[1], best):
        if p not in output:
            output[p] = {}
        output[p][k] = v
        with open(file,"w") as yaml_file:
            yaml.dump(output, yaml_file, default_flow_style=False)

def optfunc(inputargs, *args, **kwargs):
    global simplex_iteration
    cons = constraintsClass()
    if 'dir' in kwargs.keys():
        dir = kwargs['dir']
        del kwargs['dir']
        e, b, ee, be, l, g = genesisBeamFile.optfunc(inputargs, dir=dir, *args, **kwargs)
    else:
        e, b, ee, be, l, g = genesisBeamFile.optfunc(inputargs, dir='simplex/simplex_iteration_'+str(simplex_iteration), *args, **kwargs)
    simplex_iteration += 1
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
        'brightness': {'type': 'greaterthan', 'value': brightness, 'limit': 0.15, 'weight': 6},
        'bandwidth': {'type': 'lessthan', 'value': b, 'limit': 0.25, 'weight': 3},
        'pulse_energy': {'type': 'greaterthan', 'value': e, 'limit': 120, 'weight': 4},
        'bandwidth_end': {'type': 'lessthan', 'value': 1e2*abs(g.spectrum_lamdwidth_std[pos20m]), 'limit': 0.4, 'weight': 1},
        'pulse_energy_end': {'type': 'greaterthan', 'value': 1e6*abs(g.energy[pos20m]), 'limit': 200, 'weight': 2},
        'max_brightness_position': {'type': 'lessthan', 'value': abs(l), 'limit': 12, 'weight': 2.5},
        'energy': {'type': 'greaterthan', 'value': abs(momentum), 'limit': 240, 'weight': 2},
    }
    fitvalue = cons.constraints(constraintsList)
    print cons.constraintsList(constraintsList)
    print 'fitvalue[', simplex_iteration-1, '] = ', fitvalue
    saveState(inputargs, simplex_iteration-1, fitvalue)
    return fitvalue

# ******************************************************************************
POST_INJECTOR = True
CREATE_BASE_FILES = False
scaling = 5
if POST_INJECTOR and CREATE_BASE_FILES:
    optfunc = partial(optfunc, scaling=scaling, post_injector=POST_INJECTOR, verbose=False, basefiles='../basefiles_'+str(scaling)+'/')
else:
    optfunc = partial(optfunc, scaling=scaling, post_injector=POST_INJECTOR, verbose=False)

# ******************************************************************************

injector_startingvalues = [-9, 0.345, 21000000.0, -16, 0.05, -0.05]
startingvalues = best = [1.91163015e+07, -2.65517825e+01,  2.78708167e+07, -1.12304138e+00,
  1.40796310e+07,  1.66629044e+02,  3.27872006e+07,  4.22736290e+01,
 -1.48395605e-01]

injector_startingvalues = [-9, 0.345, 21000000.0, -16, 0.05, -0.05]
startingvalues = best = [2.971579612048883e7,-23.33393532244765,2.757291895432957e7,
                        -8.082040433841438,2.4141282552160002e7,183.2068993095959,
                        3.1827793513493903e7,44.76615281045321,-0.12230302685279278,0.99]

if not POST_INJECTOR:
    best = injector_startingvalues + best
elif CREATE_BASE_FILES:
    for i in [scaling]:
        pass
        # optfunc(injector_startingvalues + best, scaling=scaling, post_injector=False, verbose=False, runGenesis=False, dir='simplex/basefiles_'+str(i))

print 'best = ', best
# print 'start fitness = ', optfunc(best, dir=os.getcwd()+'/CLARA_best_simplex_elegantgenesis', scaling=5, overwrite=True, verbose=True, summary=False, post_injector=True)
# exit()

with open('simplex/best_solutions_running.csv','w') as out:
    simplex_iteration = 0
res = minimize(optfunc, best, method='nelder-mead', options={'disp': True, 'maxiter': 300, 'adaptive': True})
print res.x
#
# try:
#     print 'best fitness = ', optfunc(res.x, dir=os.getcwd()+'/simplex/final', scaling=6, overwrite=True, verbose=True, summary=True, post_injector=POST_INJECTOR)
# except:
#     pass

# saveParameterFile(res.x)
