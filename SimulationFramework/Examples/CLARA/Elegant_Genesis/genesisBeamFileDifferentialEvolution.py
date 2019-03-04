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
from scipy.optimize import differential_evolution

def saveState(args, fitness):
    with open('best_solutions_running_simplex_elegantgenesis.csv','a') as out:
        csv_out=csv.writer(out)
        args=list(args)
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
        e, b, ee, be, l, g = genesisBeamFile.optfunc(inputargs, dir='de/iteration_'+str(simplex_iteration), *args, **kwargs)
    simplex_iteration += 1
    e = 1e2*e
    ee = 1e2*ee
    pos12m = list(g.z).index(10.)
    pos20m = list(g.z).index(20.)
    constraintsList = {
        'bandwidth': {'type': 'lessthan', 'value': b, 'limit': 0.25, 'weight': 3},
        'pulse_energy': {'type': 'greaterthan', 'value': e, 'limit': 120, 'weight': 4},
        'bandwidth_end': {'type': 'lessthan', 'value': 1e2*abs(g.spectrum_lamdwidth_std[pos20m]), 'limit': 0.4, 'weight': 1},
        'pulse_energy_end': {'type': 'greaterthan', 'value': 1e6*abs(g.energy[pos20m]), 'limit': 250, 'weight': 2},
        'max_brightness_position': {'type': 'lessthan', 'value': abs(l), 'limit': 10, 'weight': 2},
    }
    fitvalue = cons.constraints(constraintsList)
    print cons.constraintsList(constraintsList)
    print 'fitvalue[', simplex_iteration-1, '] = ', fitvalue
    saveState(inputargs, fitvalue)
    return fitvalue

# ******************************************************************************
POST_INJECTOR = True
CREATE_BASE_FILES = True
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

injector_startingvalues = [-8.906156010951616,0.3420474160090586,2.0515744815221354e7,-16.281405933324855,0.05036027437405955,-0.0502414403704962]
startingvalues = best = [1.879700695321506e7,-27.832442932602543,2.583534801151637e7,-1.2644316016467563,1.446505414414888e7,181.03888866546697,
                         3.1987431329329092e7,44.128256932519484,-0.1520424081528136]

if not POST_INJECTOR:
    best = injector_startingvalues + best
elif CREATE_BASE_FILES:
    for i in [scaling]:
        pass
        # optfunc(injector_startingvalues + best, scaling=scaling, post_injector=False, verbose=False, runGenesis=False, dir='simplex/basefiles_'+str(i))

print 'best = ', best

def rangeFunc(i):
    if abs(i) > 0:
        return [0.95 * i, 1.05 * i]
    else:
        return [-1,1]

startranges = [rangeFunc(i) for i in best]

MIN = [0, -90, 0, -90, 0, 90, 0, -90, -0.2]
MAX = [33e6, 90, 33e6, 90, 45e6, 270, 32e6, 90, -0.05]

# bounds = zip(MIN,MAX)
bounds = startranges

with open('best_solutions_running_simplex_elegantgenesis.csv','w') as out:
    simplex_iteration = 0
res = differential_evolution(optfunc, bounds, maxiter=300)
print res.x
#
# try:
#     print 'best fitness = ', optfunc(res.x, dir=os.getcwd()+'/simplex/final', scaling=6, overwrite=True, verbose=True, summary=True, post_injector=POST_INJECTOR)
# except:
#     pass

# saveParameterFile(res.x)
