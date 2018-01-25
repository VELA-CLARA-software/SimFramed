from ASTRAInjector import *
import numpy as np
from constraints import *
import os
import read_beam_file as raf
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

class fitnessFunc():

    def __init__(self, args, tempdir, npart=10000, ncpu=6, overwrite=True, verbose=False):
        self.tmpdir = tempdir
        self.verbose = verbose
        linac1field, linac1phase, linac2field, linac2phase, linac3field, linac3phase, fhcfield, fhcphase, linac4field, linac4phase = args
        self.parameters = args
        self.linacfields = [linac1field, linac2field, linac3field, linac4field]
        self.dirname = os.path.basename(self.tmpdir)
        astra = ASTRAInjector(self.dirname, overwrite=overwrite)
        if not os.name == 'nt':
            astra.defineASTRACommand(['mpiexec','-np',str(ncpu),'/opt/ASTRA/astra_MPICH2.sh'])
        astra.loadSettings('short_240_12b3.settings')
        astra.modifySetting('linac1_field', abs(linac1field))
        astra.modifySetting('linac1_phase', linac1phase)
        astra.modifySetting('linac2_field', abs(linac2field))
        astra.modifySetting('linac2_phase', linac2phase)
        astra.modifySetting('linac3_field', abs(linac3field))
        astra.modifySetting('linac3_phase', linac3phase)
        astra.modifySetting('4hc_field', abs(fhcfield))
        astra.modifySetting('4hc_phase', fhcphase)
        astra.modifySetting('linac4_field', abs(linac4field))
        astra.modifySetting('linac4_phase', linac4phase)
        astra.createInitialDistribution(npart=npart, charge=250)
        astra.applySettings()
        astra.runASTRAFiles()
        # ft = feltools(self.dirname)
        # sddsfile = ft.convertToSDDS('test.5.3819.001')
        self.cons = constraintsClass()

    def calculateBeamParameters(self):
        beam = raf.beam()
        beam.slice_length = 0.1e-12
        beam.read_astra_beam_file(self.dirname+'/test.5.4936.001')
        sigmat = 1e12*np.std(beam.t)
        sigmap = np.std(beam.p)
        meanp = np.mean(beam.p)
        emitx = 1e6*beam.normalized_horizontal_emittance
        emity = 1e6*beam.normalized_vertical_emittance
        fitp = 100*sigmap/meanp
        fhcfield = self.parameters[6]
        peakI, peakIMomentumSpread, peakIEmittanceX, peakIEmittanceY, peakIMomentum = beam.sliceAnalysis()
        constraintsList = {
            'peakI': {'type': 'greaterthan', 'value': abs(peakI), 'limit': 400, 'weight': 20},
            'peakIMomentumSpread': {'type': 'lessthan', 'value': peakIMomentumSpread, 'limit': 0.1, 'weight': 3},
            'peakIEmittanceX': {'type': 'lessthan', 'value': 1e6*peakIEmittanceX, 'limit': 0.25, 'weight': 2.5},
            'peakIEmittanceY': {'type': 'lessthan', 'value': 1e6*peakIEmittanceY, 'limit': 0.25, 'weight': 2.5},
            'peakIMomentum': {'type': 'greaterthan','value': 1e-6*peakIMomentum, 'limit': 210, 'weight': 10},
            'linac fields': {'type': 'lessthan', 'value': self.linacfields, 'limit': 32, 'weight': 100},
            '4hc field': {'type': 'lessthan', 'value': fhcfield, 'limit': 35, 'weight': 100},
            'horizontal emittance': {'type': 'lessthan', 'value': emitx, 'limit': 1, 'weight': 4},
            'vertical emittance': {'type': 'lessthan', 'value': emity, 'limit': 1, 'weight': 4},
        }
        fitness = self.cons.constraints(constraintsList)
        if self.verbose:
            print self.cons.constraintsList(constraintsList)
        return fitness

def optfunc(args, dir=None, **kwargs):
    if dir == None:
        with TemporaryDirectory(dir=os.getcwd()) as tmpdir:
            fit = fitnessFunc(args, tmpdir, **kwargs)
            fitvalue = fit.calculateBeamParameters()
    else:
            fit = fitnessFunc(args, dir, **kwargs)
            fitvalue = fit.calculateBeamParameters()
    return (fitvalue,)

# if not os.name == 'nt':
#     os.chdir('/home/jkj62/ASTRAFramework/Simple')

best = [20.715615124317853,-19.490448811248754,31.499129685645897,-27.463243697074024,16.912367920177985,-12.668172100137667,27.842335415412833,188.38675572505772,20.577934659165386,3.017649550471842]
# print optfunc(best, dir=os.getcwd()+'/best_2013', npart=10000, ncpu=20, overwrite=True, verbose=True)
# exit()

startranges = [[10, 32], [-40,40], [10, 32], [-40,40], [10, 32], [-40,40], [10, 32], [135,200], [10, 32], [-40,40]]
startranges = [[0.8*i, 1.4*i] for i in best]
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

toolbox.register("evaluate", optfunc, npart=1000)
toolbox.register("mate", tools.cxBlend, alpha=0.2)
toolbox.register("mutate", tools.mutGaussian, mu=0, sigma=3, indpb=0.3)
toolbox.register("select", tools.selTournament, tournsize=3)


if __name__ == "__main__":
    random.seed(64)

    # Process Pool of 4 workers
    if not os.name == 'nt':
        pool = multiprocessing.Pool(processes=10)
    else:
        pool = multiprocessing.Pool(processes=1)
    toolbox.register("map", pool.map)

    if not os.name == 'nt':
        pop = toolbox.population(n=30)
    else:
        pop = toolbox.population(n=18)
    hof = tools.HallOfFame(10)
    stats = tools.Statistics(lambda ind: ind.fitness.values)
    stats.register("avg", numpy.mean)
    stats.register("std", numpy.std)
    stats.register("min", numpy.min)
    stats.register("max", numpy.max)

    pop, logbook = algorithms.eaSimple(pop, toolbox, cxpb=0.5, mutpb=0.2, ngen=50,
                            stats=stats, halloffame=hof)

    # print 'pop = ', pop
    print logbook
    print hof

    with open('best_solutions.csv','wb') as out:
        csv_out=csv.writer(out)
        for row in hof:
            csv_out.writerow(row)
    pool.close()

    print 'best fitness = ', optfunc(hof[0], dir=os.getcwd()+'/best', npart=50000, ncpu=40)
