import os, errno, sys
import numpy as np
import random
from scipy.constants import c
from ocelot.S2E_STFC import FEL_simulation_block
from ocelot.adaptors.genesis import generate_input, get_genesis_launcher, run_genesis, rematch_edist, edist2beam
from ocelot.gui.genesis_plot import fwhm3
from ocelot.common.math_op import *
from ocelot.cpbd.beam import Twiss
sys.path.append(os.path.abspath(__file__+'/../../../../../../'))
import SimulationFramework.Modules.read_beam_file as rbf
from SimulationFramework.Modules.optimiser import optimiser
opt = optimiser()
from SimulationFramework.Modules.constraints import constraintsClass
import matplotlib.pyplot as plt
import time
import csv
import multiprocessing
from scoop import futures
from deap import base, creator, tools, algorithms
import copy
sys.path.append('./../')
import genesisBeamFile
import id_number as idn
import id_number_server as idnserver

def saveState(args, n, params, fitness):
    with open('MOGA/best_solutions_running.csv','a') as out:
        csv_out=csv.writer(out)
        args=list(args)
        for p in params:
            args.append(p)
        args.append(n)
        args.append(fitness)
        csv_out.writerow(args)

injector_parameter_names = [
    ['CLA-HRG1-GUN-CAV', 'phase'],
    ['CLA-HRG1-GUN-SOL', 'field_amplitude'],
    ['CLA-L01-CAV', 'field_amplitude'],
    ['CLA-L01-CAV', 'phase'],
    ['CLA-L01-CAV-SOL-01', 'field_amplitude'],
    ['CLA-L01-CAV-SOL-02', 'field_amplitude'],
]
parameter_names = [
    ['CLA-L02-CAV', 'field_amplitude'],
    ['CLA-L02-CAV', 'phase'],
    ['CLA-L03-CAV', 'field_amplitude'],
    ['CLA-L03-CAV', 'phase'],
    ['CLA-L4H-CAV', 'field_amplitude'],
    ['CLA-L4H-CAV', 'phase'],
    ['CLA-L04-CAV', 'field_amplitude'],
    ['CLA-L04-CAV', 'phase'],
    ['bunch_compressor', 'set_angle'],
    ['CLA-S07-DCP-01', 'factor'],
]

POST_INJECTOR = True
startingvalues = best = np.array([31635579.78, -21.9040218, 28403832.7, -6.18551054, 24210042.9, 184.393421, 32589069.4, 44.386823, -0.124275132, 0.929080893])

def rangeFunc(i):
    if abs(i) > 0:
        return [0.95 * i, 1.05 * i]
    else:
        return [-1,1]

startranges = [rangeFunc(i) for i in best]
MIN = [0, -90, 0, -90, 0, 90, 0, -90, -0.2, 0.0]
MAX = [33e6, 90, 33e6, 90, 45e6, 270, 32e6, 90, -0.05, 3]

def checkBounds(min, max):
    def decorator(func):
        def wrapper(*args, **kargs):
            offspring = func(*args, **kargs)
            for child in offspring:
                for i in xrange(len(child)):
                    if child[i] > max[i]:
                        child[i] = max[i]
                    elif child[i] < min[i]:
                        child[i] = min[i]
            return offspring
        return wrapper
    return decorator

generateHasBeenCalled = False
def generate():
    global generateHasBeenCalled
    if not generateHasBeenCalled:
        generateHasBeenCalled = True
        return creator.Individual(list(best))
    else:
        return creator.Individual(random.uniform(a,b) for a,b in startranges)

def MOGAoptFunc(inputargs, *args, **kwargs):
    idclient = idn.zmqClient()
    n =  idclient.get_id()
    if not POST_INJECTOR:
        parameters = injector_parameter_names + parameter_names
    else:
        parameters = parameter_names
    inputlist = map(lambda a: a[0]+[a[1]], zip(parameters, inputargs))
    e, b, ee, be, l, g = genesisBeamFile.optfunc(inputlist, *args, subdir='MOGA', dir='MOGA/iteration_'+str(n), **kwargs)
    fitness = -1.0*e/b
    saveState(inputargs, n, [e, b, ee, be, l], fitness)
    return e, b, ee, be


creator.create("Fitness", base.Fitness, weights=(1.0, -1.0, 1.0, -1.0,))
creator.create("Individual", list, fitness=creator.Fitness)

toolbox = base.Toolbox()

# Attribute generator
toolbox.register("attr_bool", generate)

# Structure initializers
toolbox.register("Individual", generate)
toolbox.register("population", tools.initRepeat, list, toolbox.Individual)

if os.name == 'nt':
    toolbox.register("evaluate", MOGAoptFunc, scaling=3, post_injector=True)
else:
    toolbox.register("evaluate", MOGAoptFunc, scaling=5, post_injector=True)
toolbox.register("mate", tools.cxUniform, indpb=0.3)
toolbox.register("mutate", tools.mutGaussian, mu=0, sigma=[1e6,2,1e6,2,2e6,2,1e6,2,0.003,0.1], indpb=0.3)

toolbox.decorate("mate", checkBounds(MIN, MAX))
toolbox.decorate("mutate", checkBounds(MIN, MAX))

toolbox.register("select", tools.selNSGA2)

global_best = 0

if __name__ == "__main__":
    random.seed(43065)
    server = idnserver.zmqServer()
    server.daemon = True
    server.start()
    idclient = idn.zmqClient()
    idclient.reset_id()
    out = open('MOGA/best_solutions_running.csv','wb', buffering=0)
    # genesisBeamFile.csv_out = csv.writer(out)

    NGEN = 200
    MU = 24
    LAMBDA = 48
    CXPB = 0.7
    MUTPB = 0.2

    # Process Pool of 4 workers
    if not os.name == 'nt':
        pool = multiprocessing.Pool(processes=8)
    else:
        pool = multiprocessing.Pool(processes=2)
    toolbox.register("map", pool.map)

    if not os.name == 'nt':
        pop = toolbox.population(n=MU)
    else:
        pop = toolbox.population(n=MU)
    hof = tools.ParetoFront()
    stats = tools.Statistics(lambda ind: ind.fitness.values)
    stats.register("avg", np.mean, axis=0)
    stats.register("std", np.std, axis=0)
    stats.register("min", np.min, axis=0)
    stats.register("max", np.max, axis=0)

    opt.eaMuPlusLambda(pop, toolbox, MU, LAMBDA, CXPB, MUTPB, NGEN, stats,
                        hoffile='MOGA/CLARA_HOF_longitudinal_Genesis_DCP.csv',
                        halloffame=hof)

    pool.close()
    # print 'pop = ', pop
    # print logbook
    print hof

    # try:
    #     print 'best fitness = ', optfunc(hof[0], dir=os.getcwd()+'/CLARA_best_longitudinal_Genesis', scaling=5, overwrite=True, verbose=True, summary=True)
    #     with open('CLARA_best_longitudinal_Genesis/CLARA_best_longitudinal_Genesis_solutions.csv','wb') as out:
    #         csv_out=csv.writer(out)
    #         for row in hof:
    #             csv_out.writerow(row)
    #     with open('CLARA_best_longitudinal_Genesis/CLARA_best_longitudinal_Genesis_stats.csv','wb') as out:
    #         csv_out=csv.writer(out)
    #         for row in stats:
    #             csv_out.writerow(row)
    # except:
    with open('MOGA/CLARA_best_longitudinal_Genesis_solutions_final_DCP.csv','wb') as out:
        hof_out=csv.writer(out)
        for row in hof:
            hof_out.writerow(row)
