from ASTRAInjector import *
import numpy as np
from constraints import *
import os
import tempfile
import copy
import read_twiss_file as rtf
from deap import algorithms
from deap import base
from deap import creator
from deap import tools
import multiprocessing
from scoop import futures
import operator
import random

import csv

twiss = rtf.twiss()

astra = ASTRAInjector('.', overwrite=False)
astra.loadSettings('short_240.settings')
lines = {}
for f, s in astra.fileSettings.iteritems():
    lines[f] = astra.readFile(f)

def merge_two_dicts(x, y):
    """Given two dicts, merge them into a new dict as a shallow copy."""
    z = x.copy()
    z.update(y)
    return z

import shutil
import uuid
class TemporaryDirectory(object):
    """Context manager for tempfile.mkdtemp() so it's usable with "with" statement."""

    def __init__(self, *args, **kwargs):
        self.args = args
        self.kwargs = kwargs

    def name(self):
        return 'tmp'+str(uuid.uuid4())

    def __enter__(self, dir=os.getcwd()):
        exists = True
        while exists:
            self.name = dir + '/' + self.name()
            if not os.path.exists(self.name):
                exists=False
                os.makedirs(self.name)
        return self.name

    def __exit__(self, exc_type, exc_value, traceback):
        shutil.rmtree(self.name)

class fitnessFunc():

    def __init__(self, args, tempdir, npart=1000, ncpu=4, overwrite=True, verbose=False):
        self.tmpdir = tempdir
        self.verbose = verbose
        self.parameters = list(args)
        self.dirname = os.path.basename(self.tmpdir)
        astra = ASTRAInjector(self.dirname, overwrite=overwrite)
        if not os.name == 'nt':
            astra.defineASTRACommand(['mpiexec','-np',str(ncpu),'/opt/ASTRA/astra_MPICH2.sh'])
        else:
            astra.defineASTRACommand(['astra'])
        astra.loadSettings('short_240.settings')
        astra.fileSettings['test.2']['quad_K'] = self.parameters[0:6]
        astra.fileSettings['test.3']['quad_K'] = self.parameters[6:14]
        astra.fileSettings['test.4']['quad_K'] = self.parameters[14:16]
        astra.fileSettings['test.5']['quad_K'] = self.parameters[16:]
        # astra.fileSettings['test.2']['starting_distribution'] = astra.fileSettings['test.2']['starting_distribution']
        # astra.fileSettings['test.3']['starting_distribution'] = astra.fileSettings['test.3']['starting_distribution']
        # astra.fileSettings['test.4']['starting_distribution'] = astra.fileSettings['test.4']['starting_distribution']
        # astra.fileSettings['test.5']['starting_distribution'] = astra.fileSettings['test.5']['starting_distribution']
        # print 'Creating Initial Distribution in folder:', self.tmpdir
        astra.createInitialDistribution(npart=npart, charge=250)
        # print 'Apply Settings in folder:', self.tmpdir
        astra.applySettings(filelines=lines)
        # print 'Run ASTRA in folder:', self.tmpdir
        astra.runASTRAFiles()
        # ft = feltools(self.dirname)
        # sddsfile = ft.convertToSDDS('test.in.128.4929.128')
        # print 'Analysing Constraints in folder:', self.tmpdir
        self.cons = constraintsClass()

    def calculateTwissParameters(self):
        constraintsList = {}
        twiss.read_astra_emit_files(self.dirname+'/test.2.Zemit.001')
        constraintsList2 = {
            'max_xrms_2': {'type': 'lessthan', 'value': 1e3*max(twiss['sigma_x']), 'limit': 1, 'weight': 10},
            'max_yrms_2': {'type': 'lessthan', 'value': 1e3*max(twiss['sigma_y']), 'limit': 1, 'weight': 10},
            'min_xrms_2': {'type': 'greaterthan', 'value': 1e3*min(twiss['sigma_x']), 'limit': 0.2, 'weight': 10},
            'min_yrms_2': {'type': 'greaterthan', 'value': 1e3*min(twiss['sigma_y']), 'limit': 0.2, 'weight': 10},
            'last_exn_2': {'type': 'lessthan', 'value': 1e6*twiss['enx'][-1], 'limit': 0.8, 'weight': 100},
            'last_eyn_2': {'type': 'lessthan', 'value': 1e6*twiss['eny'][-1], 'limit': 0.8, 'weight': 100},
            'beta_x_2': {'type': 'lessthan', 'value': twiss['beta_x'], 'limit': 50, 'weight': 10},
            'beta_y_2': {'type': 'lessthan', 'value': twiss['beta_y'], 'limit': 50, 'weight': 10},
        }
        constraintsList = merge_two_dicts(constraintsList, constraintsList2)

        twiss.read_astra_emit_files(self.dirname+'/test.3.Zemit.001')
        constraintsList3 = {
            'max_xrms_3': {'type': 'lessthan', 'value': 1e3*max(twiss['sigma_x']), 'limit': 1, 'weight': 10},
            'max_yrms_3': {'type': 'lessthan', 'value': 1e3*max(twiss['sigma_y']), 'limit': 1, 'weight': 10},
            'min_xrms_3': {'type': 'greaterthan', 'value': 1e3*min(twiss['sigma_x']), 'limit': 0.2, 'weight': 10},
            'min_yrms_3': {'type': 'greaterthan', 'value': 1e3*min(twiss['sigma_y']), 'limit': 0.2, 'weight': 10},
            'last_exn_3': {'type': 'lessthan', 'value': 1e6*twiss['enx'][-1], 'limit': 0.8, 'weight': 100},
            'last_eyn_3': {'type': 'lessthan', 'value': 1e6*twiss['eny'][-1], 'limit': 0.8, 'weight': 100},
            'beta_x_3': {'type': 'lessthan', 'value': twiss['beta_x'], 'limit': 50, 'weight': 10},
            'beta_y_3': {'type': 'lessthan', 'value': twiss['beta_y'], 'limit': 50, 'weight': 10},
        }
        constraintsList = merge_two_dicts(constraintsList, constraintsList3)

        twiss.read_astra_emit_files(self.dirname+'/test.4.Zemit.001')
        constraintsList4 = {
            'min_xrms_4': {'type': 'greaterthan', 'value': 1e3*min(twiss['sigma_x']), 'limit': 0.1, 'weight': 30},
            'min_yrms_4': {'type': 'greaterthan', 'value': 1e3*min(twiss['sigma_y']), 'limit': 0.2, 'weight': 50},
            'last_yrms_4': {'type': 'greaterthan', 'value': 1e3*twiss['sigma_y'][-1], 'limit': 0.4, 'weight': 30},
            'last_xrms_4': {'type': 'lessthan', 'value': 1e3*twiss['sigma_x'][-1], 'limit': 0.2, 'weight': 30},
            'last_exn_4': {'type': 'lessthan', 'value': 1e6*twiss['enx'][-1], 'limit': 0.8, 'weight': 100},
            'last_eyn_4': {'type': 'lessthan', 'value': 1e6*twiss['eny'][-1], 'limit': 0.8, 'weight': 100},
            'beta_x_4': {'type': 'lessthan', 'value': twiss['beta_x'], 'limit': 50, 'weight': 10},
            'beta_y_4': {'type': 'lessthan', 'value': twiss['beta_y'], 'limit': 50, 'weight': 10},
        }
        constraintsList = merge_two_dicts(constraintsList, constraintsList4)

        twiss.read_astra_emit_files(self.dirname+'/test.5.Zemit.001')
        constraintsList5 = {
            'min_xrms_5': {'type': 'greaterthan', 'value': 1e3*min(twiss['sigma_x']), 'limit': 0.1, 'weight': 20},
            'min_yrms_5': {'type': 'greaterthan', 'value': 1e3*min(twiss['sigma_y']), 'limit': 0.2, 'weight': 20},
            'last_alpha_x_5': {'type': 'lessthan', 'value': abs(twiss['alpha_x'][-1]), 'limit': 2, 'weight': 10},
            'last_alpha_y_5': {'type': 'lessthan', 'value': abs(twiss['alpha_y'][-1]), 'limit': 2, 'weight': 10},
            'last_beta_x_5': {'type': 'lessthan', 'value': twiss['beta_x'], 'limit': 50, 'weight': 2},
            'tdc_phase_advance': {'type': 'equalto', 'value': twiss.interpolate(46.4378,'muy') - twiss.interpolate(40.8012,'muy'), 'limit': 0.25, 'weight': 500},
            'tdc_beta_y_greaterthan': {'type': 'greaterthan', 'value': twiss.interpolate(40.8012, 'beta_y'), 'limit': 50, 'weight': 25},
            'tdc_beta_y_lassthan': {'type': 'lessthan', 'value': twiss.interpolate(40.8012, 'beta_y'), 'limit': 100, 'weight': 25},
            'tdc_screen_beta_y': {'type': 'greaterthan', 'value': twiss.extract_values('beta_y', 40.8012, 46.4378), 'limit': 3, 'weight': 50},
            'screen_beta_x': {'type': 'equalto', 'value': twiss.interpolate(46.4378, 'beta_x'), 'limit': 5, 'weight': 50},
            'screen_beta_y': {'type': 'equalto', 'value': twiss.interpolate(46.4378, 'beta_y'), 'limit': 5, 'weight': 50},
        }
        constraintsList = merge_two_dicts(constraintsList, constraintsList5)

        fitness = self.cons.constraints(constraintsList)
        if self.verbose:
            print self.cons.constraintsList(constraintsList)
        return fitness

def optfunc(args, dir=None, **kwargs):
    if dir == None:
        # tmpdir = TemporaryDirectory(os.getcwd())
        # tmpdirname = tmpdir.name()
        with TemporaryDirectory(os.getcwd()) as tmpdirname:
            fit = fitnessFunc(args, tmpdirname, **kwargs)
            fitvalue = fit.calculateTwissParameters()
        # close(tmpdir)
    else:
        fit = fitnessFunc(args, dir, **kwargs)
        fitvalue = fit.calculateTwissParameters()
    return (fitvalue,)

# if not os.name == 'nt':
#     os.chdir('/home/jkj62/ASTRAFramework/Simple')

best = [
    1.76571, -0.334896,
    -2.38819, 3.60297, 5., -4., -3.29565, 3.38181, 1.52921, 2.22559, -2.95347, -1.71192, -5.98555, 7.72644, -5.24599
]
best = [
-1.1734910726661183,0.3619758255160046,0.26380399688214656,0.8218417571167965,0.8040239713577457,-0.9150109311909821,0.5287735079950912,0.15183492793097536,-0.5694532268795542,-0.16365543950692527,0.2500015150193916,0.45524067925448264,0.562812978218051,-0.1962066477258797,-0.5733953992449736,0.8361401644999704,0.587359403398405,-0.470367423848161,2.4140615785340254,-0.8667801563862227,-1.0663504437623548,-1.22506114725361,1.0990306191315973,0.14619286376298724,-0.0638621062250603,-0.5548413370279897,0.7534544970319805,-0.5394912911668389,-0.0407068632947221
]
# print optfunc(best, dir=os.getcwd()+'/testing', npart=10, ncpu=20, overwrite=True, verbose=True)
# exit()

best = [0 for x in best]

startranges = [[0.8*i, 1.4*i] if abs(i) > 0 else [-1,1] for i in best]

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
    toolbox.register("evaluate", optfunc, npart=100)
else:
    toolbox.register("evaluate", optfunc, npart=1000)
toolbox.register("mate", tools.cxBlend, alpha=0.2)
toolbox.register("mutate", tools.mutGaussian, mu=0, sigma=3, indpb=0.3)
toolbox.register("select", tools.selTournament, tournsize=3)


if __name__ == "__main__":
    random.seed(64)

    # Process Pool of 4 workers
    if not os.name == 'nt':
        pool = multiprocessing.Pool(processes=12)
    else:
        pool = multiprocessing.Pool(processes=3)
    toolbox.register("map", pool.map)
    # toolbox.register("map", futures.map)

    if not os.name == 'nt':
        pop = toolbox.population(n=36)
    else:
        pop = toolbox.population(n=21)
    hof = tools.HallOfFame(10)
    stats = tools.Statistics(lambda ind: ind.fitness.values)
    stats.register("avg", numpy.mean)
    stats.register("std", numpy.std)
    stats.register("min", numpy.min)
    stats.register("max", numpy.max)

    pop, logbook = algorithms.eaSimple(pop, toolbox, cxpb=0.5, mutpb=0.2, ngen=20,
                            stats=stats, halloffame=hof)

    # print 'pop = ', pop
    print logbook
    print hof

    with open('best_solutions.csv','wb') as out:
        csv_out=csv.writer(out)
        for row in hof:
            csv_out.writerow(row)
    pool.close()

    print 'best fitness = ', optfunc(hof[0], dir=os.getcwd()+'/twiss_best', npart=10000, ncpu=40)
