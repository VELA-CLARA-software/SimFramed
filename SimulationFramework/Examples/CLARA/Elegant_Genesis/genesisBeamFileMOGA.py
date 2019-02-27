import os, errno, sys
import numpy as np
import random
from scipy.constants import c
from ocelot.S2E_STFC import FEL_simulation_block
from ocelot.adaptors.genesis import generate_input, get_genesis_launcher, run_genesis, rematch_edist, edist2beam
from ocelot.gui.genesis_plot import fwhm3
from ocelot.common.math_op import *
from ocelot.cpbd.beam import Twiss
sys.path.append(os.path.abspath(__file__+'/../../../../../'))
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
import genesisBeamFile

# startingvalues = best = [30000000.0, -23, 27000000.0, -8, 24000000.0, 184, 32000000.0, 45, -0.1185]
startingvalues = best = [2.018490050471744e7,-23.04340196585034,2.934266187158792e7,-1.7771024303105327,1.7144513765057914e7,167.20031122662812,3.185636245553371e7,41.97162063476029,-0.14363986757360986,233.91419013665944]

# print 'starting values = ', best
# optfunc(best, scaling=5, post_injector=True, verbose=True)
# exit()
# print best
def rangeFunc(i):
    if abs(i) > 0:
        return [0.95 * i, 1.05 * i]
    else:
        return [-1,1]

startranges = [rangeFunc(i) for i in best]
MIN = [0, -90, 0, -90, 0, 90, 0, -90, -0.2, 100]
MAX = [33e6, 90, 33e6, 90, 45e6, 270, 32e6, 90, -0.05, 250]


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

# print generate()

creator.create("FitnessMin", base.Fitness, weights=(1.0, -1.0, 1.0, -1.0))
creator.create("Individual", list, fitness=creator.FitnessMin)

toolbox = base.Toolbox()

# Attribute generator
toolbox.register("attr_bool", generate)

# Structure initializers
toolbox.register("Individual", generate)
toolbox.register("population", tools.initRepeat, list, toolbox.Individual)

if os.name == 'nt':
    toolbox.register("evaluate", genesisBeamFile.optfunc, scaling=3, post_injector=True)
else:
    toolbox.register("evaluate", genesisBeamFile.optfunc, scaling=6, post_injector=True)
toolbox.register("mate", tools.cxUniform, indpb=0.3)
toolbox.register("mutate", tools.mutGaussian, mu=0, sigma=[1e6,2,1e6,2,2e6,2,1e6,2,0.003,5], indpb=0.3)

toolbox.decorate("mate", checkBounds(MIN, MAX))
toolbox.decorate("mutate", checkBounds(MIN, MAX))

toolbox.register("select", tools.selNSGA2)

global_best = 0

if __name__ == "__main__":
    random.seed(43065)

    out = open('best_solutions_running_simplex_elegant_genesis.csv','wb', buffering=0)
    genesisBeamFile.csv_out = csv.writer(out)

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
                        hoffile='CLARA_HOF_longitudinal_Genesis.csv',
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
    with open('CLARA_best_longitudinal_Genesis_solutions_final.csv','wb') as out:
        hof_out=csv.writer(out)
        for row in hof:
            hof_out.writerow(row)
