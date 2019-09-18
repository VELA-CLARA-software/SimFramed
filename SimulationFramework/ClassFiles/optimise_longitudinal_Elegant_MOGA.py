import os, errno, sys
import numpy as np
import random
sys.path.append('./../../')
from SimulationFramework.Modules.optimiser import optimiser
from SimulationFramework.ClassFiles.Optimise_longitudinal_Elegant import Optimise_Elegant
opt = optimiser()
import multiprocessing
# from scoop import futures
import deap#, deap.base, deap.creator, deap.tools
from copy import copy
import csv
import SimulationFramework.Modules.id_number as idn
import SimulationFramework.Modules.id_number_server as idnserver

class MOGA(Optimise_Elegant):

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

    def __init__(self):
        super(MOGA, self).__init__()
        self.global_best = 0
        self.POST_INJECTOR = True

    def create_weights_function(self, weights=(-1.0, 1.0, -1.0, 1.0, )):
        deap.creator.create("Fitness", deap.base.Fitness, weights=weights)
        deap.creator.create("Individual", list, fitness=deap.creator.Fitness)

    def create_toolbox(self):
        self.toolbox = deap.base.Toolbox()
        # Attribute generator
        self.toolbox.register("attr_bool", self.generate)
        # Structure initializers
        self.toolbox.register("Individual", self.generate)
        self.toolbox.register("population", deap.tools.initRepeat, list, self.toolbox.Individual)

    def create_fitness_function(self, function, **kwargs):
        self.toolbox.register("evaluate", function, **kwargs)#scaling=3, post_injector=True)

    def create_mating_function(self, method, **kwargs):
        self.toolbox.register("mate", method, **kwargs)

    def create_uniform_mating_function(self, probability=0.3):
        self.create_mating_function(deap.tools.cxUniform, indpb=probability)

    def create_mutation_function(self, method, **kwargs):
        self.toolbox.register("mutate", method, **kwargs)

    def create_gaussian_mutation_function(self, probability=0.3, mu=0, sigma=[1e6,2,1e6,2,2e6,2,1e6,2,0.003,0.1]):
        self.create_mutation_function(deap.tools.mutGaussian, mu=mu, sigma=sigma, indpb=probability)

    def add_bounds(self, MIN, MAX):
        self.toolbox.decorate("mate", self.checkBounds(MIN, MAX))
        self.toolbox.decorate("mutate", self.checkBounds(MIN, MAX))

    def create_selection_function(self, method, **kwargs):
        self.toolbox.register("select", method, **kwargs)

    def create_NSGA2_selection_function(self, **kwargs):
        self.create_selection_function(deap.tools.selNSGA2, **kwargs)

    def saveState(self, args, n, params, fitness):
        with open('MOGA/best_solutions_running.csv','a') as out:
            csv_out=csv.writer(out)
            args=list(args)
            for p in params:
                args.append(p)
            args.append(n)
            args.append(fitness)
            csv_out.writerow(args)

    def rangeFunc(self, i):
        if abs(i) > 0:
            return [0.95 * i, 1.05 * i]
        else:
            return [-1,1]

    def checkBounds(self, min, max):
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

    def generate(self):
        if not self.generateHasBeenCalled:
            self.generateHasBeenCalled = True
            return deap.creator.Individual(list(self.best))
        else:
            return deap.creator.Individual(random.uniform(a,b) for a,b in self.startranges)

    def MOGAoptFunc(self, inputargs, *args, **kwargs):
        e, b, ee, be, l, g = self.OptimisingFunction(inputargs, **kwargs)
        fitness = -1.0*e/b
        print ('fitvalue[', self.opt_iteration, '] - E=', 1e2*e, '  BW=', b, '  fitness=', fitness)
        self.saveState(inputargs, self.opt_iteration, [e, b, ee, be, l], fitness)
        return e, b, e/b

    def OptimisingFunction(self, inputargs, **kwargs):
        self.optdir = 'MOGA/iteration_'
        if not self.post_injector:
            parameternames = self.injector_parameter_names + self.parameter_names
        else:
            parameternames = copy(self.parameter_names)

        self.inputlist = list(map(lambda a: a[0]+[a[1]], zip(parameternames, inputargs)))

        self.linac_fields = np.array([i[2] for i in self.inputlist if i[1] == 'field_amplitude'])
        self.linac_phases = np.array([i[2] for i in self.inputlist if i[1] == 'phase'])

        idclient = idn.zmqClient()
        n =  idclient.get_id()
        # print('id n = ', n)
        self.opt_iteration = n

        dir = self.optdir+str(self.opt_iteration)
        if not os.path.exists(dir):
            os.makedirs(dir)
        self.setup_lattice(self.inputlist, dir)
        self.before_tracking()
        fitvalue = self.track()

        if isinstance(self.opt_iteration, int):
            self.opt_iteration += 1

        if e < 0.01:
            print ('e too low! ', e)
            l = 500
        self.resultsDict.update({'e': e, 'b': b, 'ee': ee, 'be': be, 'l': l, 'g': g, 'brightness': (1e-4*e)/(1e-2*b), 'momentum': 1e-6*np.mean(g.momentum)})
        return e, b, ee, be, l, g

    def initialise_population(self, best, npop, sigma=None):
        self.best = best
        self.generateHasBeenCalled = False
        if sigma is None:
            self.startranges = [self.rangeFunc(i) for i in best]
        else:
            self.startranges = [[b-s,b+s] for b, s in zip(best,sigma)]
        self.pop = self.toolbox.population(n=npop)

    def initialise_MOGA(self, seed=6546841):
        random.seed(seed)
        self.hof = deap.tools.ParetoFront()
        self.stats = deap.tools.Statistics(lambda ind: ind.fitness.values)
        self.stats.register("avg", np.mean, axis=0)
        self.stats.register("std", np.std, axis=0)
        self.stats.register("min", np.min, axis=0)
        self.stats.register("max", np.max, axis=0)

        self.server = idnserver.zmqServer()
        self.server.daemon = True
        self.server.start()

    def eaMuPlusLambda(self, nSelect, nChildren, crossoverprobability, mutationprobability, ngenerations):
        self.optdir = 'MOGA/iteration_'
        self.best_changes = './MOGA_best_changes.yaml'
        self.bestfit = 1e26
        out = open('MOGA/best_solutions_running.csv','wb', buffering=0)
        self.csv_out = csv.writer(out)
        opt.eaMuPlusLambda(self.pop, self.toolbox, nSelect, nChildren, crossoverprobability, mutationprobability, ngenerations, self.stats,
                            hoffile='MOGA/CLARA_HOF_longitudinal_Genesis_DCP.csv', halloffame=self.hof)
