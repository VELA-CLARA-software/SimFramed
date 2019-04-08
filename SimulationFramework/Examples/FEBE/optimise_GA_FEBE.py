import os, sys
sys.path.append('../../../')
from SimulationFramework.Framework import *
import numpy as np
from SimulationFramework.Modules.constraints import *
import SimulationFramework.Modules.read_twiss_file as rtf
import SimulationFramework.Modules.read_beam_file as rbf
from SimulationFramework.Modules.optimiser import optimiser
opt = optimiser()
from deap import algorithms
from deap import base
from deap import creator
from deap import tools
import operator
import random
import multiprocessing
from scoop import futures
import csv

import shutil
import uuid

import signal
signal.signal(signal.SIGINT, opt.finish_running)

def merge_two_dicts(x, y):
    """Given two dicts, merge them into a new dict as a shallow copy."""
    z = x.copy()
    z.update(y)
    return z

class TemporaryDirectory(object):
    """Context manager for tempfile.mkdtemp() so it's usable with "with" statement."""

    def __init__(self, *args, **kwargs):
        self.args = args
        self.kwargs = kwargs

    def tempname(self):
        return 'tmp'+str(uuid.uuid4())

    def __enter__(self, dir=os.getcwd()):
        exists = True
        while exists:
            self.name = dir + '/' + self.tempname()
            if not os.path.exists(self.name):
                exists=False
                os.makedirs(self.name)
        return self.name

    def __exit__(self, exc_type, exc_value, traceback):
        shutil.rmtree(self.name)

def create_base_files(scaling):
    framework = Framework('basefiles_FEBE_'+str(scaling), overwrite=False)
    framework.loadSettings('Lattices/clara400_v12_FEBE.def')
    if not os.name == 'nt':
        framework.defineASTRACommand(['mpiexec','-np',str(3*scaling),'/opt/ASTRA/astra_MPICH2.sh'])
        framework.defineGeneratorCommand(['/opt/ASTRA/generator.sh'])
        framework.defineCSRTrackCommand(['/opt/OpenMPI-1.4.3/bin/mpiexec','-n',str(3*scaling),'/opt/CSRTrack/csrtrack_openmpi.sh'])
    framework.defineElegantCommand(['elegant'])
    framework['S07'].file_block['input']['prefix'] = '../../basefiles_'+str(scaling)+'/'
    framework.track(track=True, startfile='S07')

# for i in [6]:
#     create_base_files(i)
# exit()

framework = Framework('basefiles_6', overwrite=False)
framework.loadSettings('Lattices/clara400_v12_FEBE.def')
# posstart = 60
# parameters = [a['k1'] for a in framework.getElementType('quadrupole') if a['position_start'][2] > posstart]
variables = [a for a in framework.getElementType('quadrupole') if 'S07F' in a.objectName]
preparameters = parameters = [a['k1'] for a in variables]
print 'parameters = ', parameters
# exit()
best = parameters

class fitnessFunc():

    def __init__(self, args, tempdir, scaling=4, overwrite=True, verbose=False, summary=False, clean=False):
        global preparameters, variables
        self.cons = constraintsClass()
        self.beam = rbf.beam()
        self.twiss = rtf.twiss()
        self.tmpdir = tempdir
        self.scaling = scaling
        self.verbose = verbose
        self.summary = summary
        self.parameters = preparameters + list(args)
        self.dirname = os.path.basename(self.tmpdir)
        self.framework = Framework(self.dirname, clean=clean, verbose=False)
        self.framework.loadSettings('Lattices/clara400_v12_FEBE.def')
        if not os.name == 'nt':
            self.framework.defineASTRACommand(['mpiexec','-np',str(3*scaling),'/opt/ASTRA/astra_MPICH2.sh'])
            self.framework.defineGeneratorCommand(['/opt/ASTRA/generator.sh'])
            self.framework.defineCSRTrackCommand(['/opt/OpenMPI-1.4.3/bin/mpiexec','-n',str(3*scaling),'/opt/CSRTrack/csrtrack_openmpi.sh'])
            # self.framework.defineElegantCommand(['mpiexec','-np',str(3*scaling),'Pelegant'])
        # else:
        self.framework.defineElegantCommand(['elegant'])
        # self.framework.setElementType('quadrupole','k1', self.parameters)
        [setattr(a, 'k1', 2*a['k1']) for a in variables]


    def between(self, value, minvalue, maxvalue, absolute=True):
        if absolute:
            result = max([minvalue,min([maxvalue,abs(value)])])
        else:
            result = np.sign(value)*max([minvalue,min([maxvalue,abs(value)])])
        return result

    def calculateBeamParameters(self):
            twiss = self.twiss
        # try:
            startS = self.framework['S07'].startObject['position_start'][2]
            self.framework['S07'].file_block['input']['prefix'] = '../../basefiles_'+str(self.scaling)+'/'
            self.framework.track(track=True, startfile='S07')

            constraintsList = {}
            constraintsListQuads = {
                'max_k': {'type': 'lessthan', 'value': [abs(p) for p in self.parameters], 'limit': 2.5, 'weight': 10},

            }
            # constraintsList = merge_two_dicts(constraintsList, constraintsListQuads)
            twiss.reset_dicts()
            twiss.read_sdds_file( self.dirname+'/'+'FEBE.twi' )
            twiss.read_sdds_file( self.dirname+'/'+'FEBE.sig' )
            ip_position = self.framework['FEBE'].findS('CLA-FEB-W-FOCUS-01')[0][1]
            constraintsListS07 = {
                'dechirper_sigma_x': {'type': 'lessthan', 'value': 1e3*twiss.interpolate(ip_position, 'betax', index='s'), 'limit': 0.1, 'weight': 10},
                'dechirper_sigma_y': {'type': 'lessthan', 'value': 1e3*twiss.interpolate(ip_position, 'betay', index='s'), 'limit': 0.1, 'weight': 10},
            }
            constraintsList = merge_two_dicts(constraintsList, constraintsListS07)
            fitness = self.cons.constraints(constraintsList)
            if self.verbose:
                print self.cons.constraintsList(constraintsList)
            if self.summary:
                self.framework.createHDF5Summary(reference='Transverse_GA')
            return fitness
        # except:
        #     return 1e6

def optfunc(args, dir=None, **kwargs):
    if dir == None:
        with TemporaryDirectory(dir=os.getcwd()) as tmpdir:
            fit = fitnessFunc(args, tmpdir, **kwargs)
            fitvalue = fit.calculateBeamParameters()
    else:
            fit = fitnessFunc(args, dir, **kwargs)
            fitvalue = fit.calculateBeamParameters()
    return (fitvalue,)

print 'starting values = ', best
print optfunc(best, dir=os.getcwd()+'/transverse_FEBE_best_Short_240', scaling=6, overwrite=True, verbose=True, summary=False)
# fit = fitnessFunc(best, os.getcwd()+'/test', scaling=4, overwrite=True, verbose=True, summary=False)
# fitvalue = fit.calculateBeamParameters()
exit()


# startranges = [[10, 32], [-40,40], [10, 32], [-40,40], [10, 32], [-40,40], [10, 32], [135,200], [10, 32], [-40,40], [0.8,0.15]]
startranges = [[0.8*i, 1.2*i] if abs(i) > 0 else [-0.1,0.1] for i in best]
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
    toolbox.register("evaluate", optfunc, scaling=6)
else:
    toolbox.register("evaluate", optfunc, scaling=6)
toolbox.register("mate", tools.cxBlend, alpha=0.2)
toolbox.register("mutate", tools.mutGaussian, mu=0, sigma=3, indpb=0.3)
toolbox.register("select", tools.selTournament, tournsize=3)


if __name__ == "__main__":
    random.seed(64)

    # Process Pool of 4 workers
    if not os.name == 'nt':
        pool = multiprocessing.Pool(processes=24)
    else:
        pool = multiprocessing.Pool(processes=6)
    toolbox.register("map", pool.map)

    if not os.name == 'nt':
        pop = toolbox.population(n=24)
    else:
        pop = toolbox.population(n=24)
    hof = tools.HallOfFame(10)
    stats = tools.Statistics(lambda ind: ind.fitness.values)
    stats.register("avg", np.mean)
    stats.register("std", np.std)
    stats.register("min", np.min)
    stats.register("max", np.max)

    pop, logbook = opt.eaSimple(pop, toolbox, cxpb=0.5, mutpb=0.2, ngen=50,
                            stats=stats, halloffame=hof)

    # print 'pop = ', pop
    print logbook
    print hof

    try:
        print 'best fitness = ', optfunc(hof[0], dir=os.getcwd()+'/transverse_DCP_best_Short_240', scaling=6, overwrite=True, verbose=True, summary=False, clean=True)
        with open('transverse_DCP_best_Short_240/transverse_DCP_best_solutions.csv','wb') as out:
            csv_out=csv.writer(out)
            for row in hof:
                csv_out.writerow(row)
    except:
        with open('transverse_DCP_best_Short_240_solutions.csv.tmp','wb') as out:
            csv_out=csv.writer(out)
            for row in hof:
                csv_out.writerow(row)
    pool.close()
