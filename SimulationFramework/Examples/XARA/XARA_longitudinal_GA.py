import sys, os
from copy import copy
import numpy as np
sys.path.append('./../../../')
# from SimulationFramework.ClassFiles.Optimise_longitudinal_Elegant import Optimise_Elegant, saveState
from XARA_nelder_mead import XARA
from SimulationFramework.Modules.merge_two_dicts import merge_two_dicts
from functools import partial
from ruamel import yaml
from shutil import copyfile
import SimulationFramework.Framework as fw
from deap import algorithms
from deap import base
from deap import creator
from deap import tools
import operator
import random
import multiprocessing
from scoop import futures
import csv
from SimulationFramework.Modules.optimiser import optimiser
eaopt = optimiser()
# sys.path.append('./../CLARA/Elegant_Genesis/')
import SimulationFramework.Modules.id_number as idn
import SimulationFramework.Modules.id_number_server as idnserver
import best_value as bestclient
import best_value_server as bestserver

# import matplotlib.pyplot as plt
# plt.ion()
# plt.show()
# fig, ax = plt.subplots(2, 3, sharey=False,
#                        figsize=(16, 6))

class XARA_GA(XARA):

    def __init__(self):
        super(XARA_GA, self).__init__()

    def OptimisingFunction2(self, inputargs, *args, **kwargs):
        idclient = idn.zmqClient()
        n =  idclient.get_id()
        self.opt_iteration = n
        bestcli = bestclient.zmqClient()
        if not self.POST_INJECTOR:
            parameternames = self.injector_parameter_names + self.parameter_names
        else:
            parameternames = copy(self.parameter_names)
        self.inputlist = [a[0]+[a[1]] for a in zip(parameternames, inputargs)]

        self.linac_fields = np.array([i[2] for i in self.inputlist if i[1] == 'field_amplitude'])
        self.linac_phases = np.array([i[2] for i in self.inputlist if i[1] == 'phase'])

        if 'dir' in list(kwargs.keys()):
            dir = kwargs['dir']
            del kwargs['dir']
        else:
            dir = self.optdir+str(self.opt_iteration)

        self.setup_lattice(self.inputlist, dir, changes=self.changes, *args, **kwargs)
        self.before_tracking()
        fitvalue = self.track()
        constraintsList = self.calculate_constraints()
        fitvalue = self.cons.constraints(constraintsList)
        print('fitvalue[', self.opt_iteration, '] = ', fitvalue)
        # saveState(self.subdir, inputargs, self.opt_iteration, fitvalue)
        # print('bestcli = ', bestcli.set_best(fitvalue, self.opt_iteration))
        if fitvalue < bestcli.set_best(fitvalue, self.opt_iteration):
            print(self.cons.constraintsList(constraintsList))
            print('!!!!!!  New best = ', fitvalue, inputargs)
            copyfile(dir+'/changes.yaml', self.best_changes)
            self.bestfit = fitvalue
        # except:
        #     fitvalue = 1e12
        return (fitvalue,)

opt = XARA_GA()
opt.set_changes_file(['nelder_mead_best_changes.yaml', './transverse_best_changes.yaml'])
opt.set_lattice_file('Lattices/claraX400_v12_80MVm_Elegant.def')
opt.set_start_file('CLARAX')
best = opt.load_best('GA_best_changes.yaml')

opt.subdir = 'longitudinal_GA'
opt.optdir = opt.subdir + '/iteration_'
opt.best_changes = './GA_best_changes.yaml'

def optfunc(inputargs, **kwargs):
    try:
        return opt.OptimisingFunction2(inputargs, **kwargs)
    except Exception as e:
        print 'Error! = ', e
        return (1e22,)


# startranges = [[0.95*i, 1.05*i] if abs(i) > 0 else [-20,20] for i in best]
startranges = [[b-s, b+s] for b,s in zip(best,[5e6, 5, 5e6, 5, 5e6, 5, 5e6, 5, 5e6, 5,  5e6, 5,  5e6, 5, 0.005, 0.1, 5e-5])]
# print(('startranges = ', startranges))
generateHasBeenCalled = False
def generate():
    global generateHasBeenCalled
    if not generateHasBeenCalled:
        generateHasBeenCalled = True
        return creator.Individual(list(best))
    else:
        return creator.Individual(random.uniform(a,b) for a,b in startranges)

# print generate()

creator.create("FitnessMin", base.Fitness, weights=(-1.0,))
creator.create("Individual", list, fitness=creator.FitnessMin)

toolbox = base.Toolbox()

# Attribute generator
toolbox.register("attr_bool", generate)

# Structure initializers
toolbox.register("Individual", generate)
toolbox.register("population", tools.initRepeat, list, toolbox.Individual)

if os.name == 'nt':
    toolbox.register("evaluate", optfunc)
else:
    toolbox.register("evaluate", optfunc)

toolbox.register("mate", tools.cxBlend, alpha=0.2)
toolbox.register("mutate", tools.mutGaussian, mu=0, sigma=[5e6, 5, 5e6, 5, 5e6, 5, 5e6, 5, 5e6, 5,  5e6, 5,  5e6, 5, 0.005, 0.1, 5e-5], indpb=0.3)
toolbox.register("select", tools.selTournament, tournsize=3)


if __name__ == "__main__":

    server = idnserver.zmqServer()
    server.daemon = True
    server.start()

    bestsrv = bestserver.zmqServer()
    bestsrv.daemon = True
    bestsrv.start()

    global hof
    random.seed(64)

    out = open('best_solutions_longitudinal_GA_elegant.csv','w')
    hoffile='XARA_HOF_longitudinal_Elegant.csv'
    # FF.csv_out = csv.writer(out)

    hof = tools.HallOfFame(10)
    stats = tools.Statistics(lambda ind: ind.fitness.values)
    stats.register("avg", np.mean)
    stats.register("std", np.std)
    stats.register("min", np.min)
    stats.register("max", np.max)

    nGenerations = 200
    populationSize = 20
    nSelect = populationSize
    nChildren = 2*nSelect
    crossoverprobability = 0.6
    mutationprobability = 0.2
    ngenerations = 200

    # Process Pool of 4 workers
    if not os.name == 'nt':
        pool = multiprocessing.Pool(processes=populationSize)
    else:
        pool = multiprocessing.Pool(processes=populationSize/4)
    toolbox.register("map", pool.map)

    if not os.name == 'nt':
        pop = toolbox.population(n=populationSize)
    else:
        pop = toolbox.population(n=populationSize)


    pop, logbook = eaopt.eaMuPlusLambda(pop, toolbox, mu=nSelect, lambda_=nChildren, cxpb=crossoverprobability, mutpb=mutationprobability, ngen=ngenerations,
                            stats=stats, halloffame=hof, verbose=True,
                            hoffile=hoffile)

    pool.close()
    # print 'pop = ', pop
    print(logbook)
    print(hof)

    # try:
    print('best fitness = ', optfunc(hof[0], dir=os.getcwd()+'/XARA_best_longitudinal_elegant', scaling=6, overwrite=True, verbose=True, summary=True, post_injector=post_inj))
