import os, errno, sys
import numpy as np
import random
sys.path.append(os.path.abspath(__file__+'/../../../../../'))
from SimulationFramework.Modules.constraints import constraintsClass
from SimulationFramework.Modules.nelder_mead import nelder_mead
import time
import csv
import copy
import genesisBeamFile
from functools import partial
# from scipy.optimize import minimize
from shutil import copyfile

def saveState(args, n, fitness):
    with open('nelder_mead/best_solutions_running.csv','a') as out:
        csv_out=csv.writer(out)
        args=list(args)
        args.append(n)
        args.append(fitness)
        csv_out.writerow(args)

def saveParameterFile(best, file='clara_elegantgenesis_best.yaml'):
    if POST_INJECTOR:
        allparams = zip(*(parameternames))
    else:
        allparams = zip(*(injparameternames+parameternames))
    output = {}
    for p, k, v in zip(allparams[0], allparams[1], best):
        if p not in output:
            output[p] = {}
        output[p][k] = v
        with open(file,"w") as yaml_file:
            yaml.dump(output, yaml_file, default_flow_style=False)

def optfunc(inputargs, *args, **kwargs):
    global nelder_mead_iteration, bestfit
    cons = constraintsClass()
    if 'dir' in kwargs.keys():
        dir = kwargs['dir']
        del kwargs['dir']
    else:
        dir = 'nelder_mead/nelder_mead_iteration_'+str(nelder_mead_iteration)
    e, b, ee, be, l, g = genesisBeamFile.optfunc(inputargs, dir=dir, *args, **kwargs)
    nelder_mead_iteration += 1
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
        'bandwidth': {'type': 'lessthan', 'value': b, 'limit': 0.1, 'weight': 3},
        'pulse_energy': {'type': 'greaterthan', 'value': e, 'limit': 120, 'weight': 4},
        'bandwidth_end': {'type': 'lessthan', 'value': 1e2*abs(g.spectrum_lamdwidth_std[pos20m]), 'limit': 0.35, 'weight': 1},
        'pulse_energy_end': {'type': 'greaterthan', 'value': 1e6*abs(g.energy[pos20m]), 'limit': 250, 'weight': 2},
        'max_brightness_position': {'type': 'lessthan', 'value': abs(l), 'limit': 11, 'weight': 2.5},
        'min_brightness_position': {'type': 'greaterthan', 'value': abs(l), 'limit': 10, 'weight': 2.5},
        'energy': {'type': 'greaterthan', 'value': abs(momentum), 'limit': 240, 'weight': 2},
        'energy': {'type': 'lessthan', 'value': abs(momentum), 'limit': 250, 'weight': 2},
    }
    fitvalue = cons.constraints(constraintsList)
    print cons.constraintsList(constraintsList)
    print 'fitvalue[', nelder_mead_iteration-1, '] = ', fitvalue
    saveState(inputargs, nelder_mead_iteration-1, fitvalue)
    if fitvalue < bestfit:
        print '!!!!!!  New best = ', fitvalue
        copyfile(dir+'/changes.yaml','./nelder_mead_best_changes.yaml')
        bestfit = fitvalue
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

injector_startingvalues = [-9.,0.345,2.1e7,-16.,0.052500000000000005,-0.05]
startingvalues = best = np.array([3.003162576456073e7,-22.688094300982396,2.7858659446474683e7,-8.171745448661719,
                        2.442011657057434e7,181.8366939201224,3.2198182091182925e7,45.36727642193955,-0.12145534694213986, 1.0015993470346622])

if not POST_INJECTOR:
    best = injector_startingvalues + best
elif CREATE_BASE_FILES:
    for i in [scaling]:
        pass
        # optfunc(injector_startingvalues + best, scaling=scaling, post_injector=False, verbose=False, runGenesis=False, dir='nelder_mead/basefiles_'+str(i))

print 'best = ', best
# print 'start fitness = ', optfunc(best, dir=os.getcwd()+'/CLARA_best_nelder_mead_elegantgenesis', scaling=5, overwrite=True, verbose=True, summary=False, post_injector=True)
# exit()
bestfit = 1e26

with open('nelder_mead/best_solutions_running.csv','w') as out:
    nelder_mead_iteration = 0
res = nelder_mead(optfunc, best, step=[1e6, 5, 1e6, 5, 1e6, 5, 1e6, 5, 0.01, 0.25], max_iter=300)
print res
#
# try:
#     print 'best fitness = ', optfunc(res.x, dir=os.getcwd()+'/nelder_mead/final', scaling=6, overwrite=True, verbose=True, summary=True, post_injector=POST_INJECTOR)
# except:
#     pass

# saveParameterFile(res)
