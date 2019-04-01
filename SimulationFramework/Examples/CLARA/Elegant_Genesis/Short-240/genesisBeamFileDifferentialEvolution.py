import os, errno, sys, copy
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
import uuid
import id_number as idn
import id_number_server as idnserver

def saveState(args, n, fitness):
    with open('de/best_solutions_running.csv','a') as out:
        csv_out=csv.writer(out)
        args=list(args)
        args.append(n)
        args.append(fitness)
        csv_out.writerow(args)

def optfunc(inputargs, *args, **kwargs):
    idclient = idn.zmqClient()
    n =  idclient.get_id()#str(uuid.uuid4())
    # kwargs = args[0]
    # print 'optfunc kwargs = ', kwargs


    cons = constraintsClass()
    if 'dir' in kwargs.keys():
        dir = kwargs['dir']
        del kwargs['dir']
        e, b, ee, be, l, g = genesisBeamFile.optfunc(inputargs, dir=dir, **kwargs)
    else:
        e, b, ee, be, l, g = genesisBeamFile.optfunc(inputargs, dir='de/de_iteration_'+str(n), **kwargs)
    e = 1e2*e
    ee = 1e2*ee
    pos12m = list(g.z).index(10.)
    pos20m = list(g.z).index(20.)
    momentum = 1e-6*np.mean(g.momentum)
    print 'Momentum = ', momentum
    if e < 1:
        l = 500
    constraintsList = {
        'bandwidth': {'type': 'lessthan', 'value': b, 'limit': 0.25, 'weight': 3},
        'pulse_energy': {'type': 'greaterthan', 'value': e, 'limit': 120, 'weight': 4},
        'bandwidth_end': {'type': 'lessthan', 'value': 1e2*abs(g.spectrum_lamdwidth_std[pos20m]), 'limit': 0.4, 'weight': 1},
        'pulse_energy_end': {'type': 'greaterthan', 'value': 1e6*abs(g.energy[pos20m]), 'limit': 250, 'weight': 2},
        'max_brightness_position': {'type': 'lessthan', 'value': abs(l), 'limit': 10, 'weight': 2},
        'energy': {'type': 'greaterthan', 'value': abs(momentum), 'limit': 240, 'weight': 2},
    }
    fitvalue = cons.constraints(constraintsList)
    print cons.constraintsList(constraintsList)
    print 'fitvalue[', n, '] = ', fitvalue
    try:
        saveState(inputargs, n, fitvalue)
    except Exception as e:
        print 'Error! ', e
    return fitvalue


if __name__ == '__main__':

    print 'in main!'
    server = idnserver.zmqServer()
    server.daemon = True
    server.start()
    idclient = idn.zmqClient()
    idclient.reset_id()
    print 'id reset'
    # ******************************************************************************
    POST_INJECTOR = True
    CREATE_BASE_FILES = True
    scaling = 5
    # if POST_INJECTOR and CREATE_BASE_FILES:
    #     optfunc = partial(optfunc, scaling=scaling, post_injector=POST_INJECTOR, verbose=False, basefiles='../basefiles_'+str(scaling)+'/')
    # else:
    #     optfunc = partial(optfunc, scaling=scaling, post_injector=POST_INJECTOR, verbose=False)

    # ******************************************************************************

    injector_startingvalues = [-9, 0.345, 21000000.0, -16, 0.05, -0.05]
    startingvalues = best = [18098078.513507977,-28.111766926229137,31441717.849741504,-1.0448097057137171,14144444.379584715,
                             168.92627603174682,31641201.21981612,45.08581803817373,-0.1570956730945702]

    # injector_startingvalues = [-8.906156010951616,0.3420474160090586,2.0515744815221354e7,-16.281405933324855,0.05036027437405955,-0.0502414403704962]
    # startingvalues = best = [19026408.665955413,-27.402785729222668,24380076.842652954,-1.2693696634519327,15884355.067580678,158.63744077479666,
                             # 29046200.15473125,41.711222422668044,-0.15134226422592262]

    if not POST_INJECTOR:
        best = injector_startingvalues + best
    elif CREATE_BASE_FILES:
        for i in [scaling]:
            pass
            # optfunc(injector_startingvalues + best, scaling=scaling, post_injector=False, verbose=False, runGenesis=False, dir='de/basefiles_'+str(i))

    print 'best = ', best

    def rangeFunc(n, i):
        if abs(i) > 0:
            ans = [0.9 * i, 1.1 * i]
            ans[0] = ans[0] if ans[0] > MIN[n] else MIN[n]
            ans[1] = ans[1] if ans[1] < MAX[n] else MAX[n]
            return ans
        else:
            return [MIN[n], MAX[n]]

    MIN = [0, -90, 0, -90, 0, 90, 0, -90, -0.2]
    MAX = [33e6, 90, 33e6, 90, 45e6, 270, 32e6, 90, -0.05]

    startranges = [rangeFunc(n, i) for n, i in enumerate(best)]

    # bounds = zip(MIN,MAX)
    bounds = startranges
    print 'bounds = ', bounds

    # if using basefiles
    # , 'basefiles': '../basefiles_'+str(scaling)+'/'

    if os.path.isfile('de/best_solutions_running.csv'):
        os.remove('de/best_solutions_running.csv')
    with open('de/best_solutions_running.csv','w') as out:
        simplex_iteration = 0
    res = differential_evolution(optfunc, bounds, args={'scaling': scaling, 'post_injector': POST_INJECTOR, 'verbose': False},
                                 maxiter=300, workers=8, updating='deferred')
    print res.x
    #
    # try:
    #     print 'best fitness = ', optfunc(res.x, dir=os.getcwd()+'/simplex/final', scaling=6, overwrite=True, verbose=True, summary=True, post_injector=POST_INJECTOR)
    # except:
    #     pass

    # saveParameterFile(res.x)
