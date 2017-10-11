from ASTRAInjector import *
import blackbox as bb
# import sdds
import numpy as np
from constraints import *
import tempfile
import os
import shutil
import read_astra_file as raf
from deap import base
from deap import benchmarks
from deap import creator
from deap import tools
import operator
import random

class fitnessFunc():

    def __init__(self, args):
        linac1field, linac1phase, linac2field, linac2phase, linac3field, linac3phase, fhcfield, fhcphase, linac4field, linac4phase = args
        # print 'args = ', args
        self.parameters = args
        self.linacfields = [linac1field, linac2field, linac3field, linac4field]
        self.tmpdir = tempfile.mkdtemp(dir=os.getcwd())
        self.dirname = os.path.basename(self.tmpdir)
        astra = ASTRAInjector(self.dirname, overwrite=True)
        if not os.name == 'nt':
            astra.defineASTRACommand(['mpiexec','-np','12','/opt/ASTRA/astra_MPICH2_Quiet.sh'])
        astra.loadSettings('short_240.settings')
        astra.modifySetting('linac1_field', linac1field)
        astra.modifySetting('linac1_phase', linac1phase)
        astra.modifySetting('linac2_field', linac2field)
        astra.modifySetting('linac2_phase', linac2phase)
        astra.modifySetting('linac3_field', linac3field)
        astra.modifySetting('linac3_phase', linac3phase)
        astra.modifySetting('4hc_field', fhcfield)
        astra.modifySetting('4hc_phase', fhcphase)
        astra.modifySetting('linac4_field', linac4field)
        astra.modifySetting('linac4_phase', linac4phase)
        astra.createInitialDistribution(npart=100, charge=250)
        astra.applySettings()
        astra.runASTRAFiles()
        ft = feltools(self.dirname)
        sddsfile = ft.convertToSDDS('test.in.128.4929.128')

    def removeTempDirectory(self):
        shutil.rmtree(self.tmpdir)

    def calculateBeamParameters(self):
        beam = raf.read_astra_beam_file(self.dirname+'/test.in.128.4929.128')
        c = [0] * 5
        w = [2,1,2,5,5]
        sigmat = 1e12*np.std(beam['t'])
        sigmap = np.std(beam['p'])
        meanp = np.mean(beam['p'])
        fitp = 100*sigmap/meanp
        fhcfield = self.parameters[6]
        c[0] = 0 if sigmat < 0.5 else (np.abs(sigmat-0.5))
        c[1] = 0 if fitp < 0.5 else (np.abs(fitp-0.5))
        c[2] = 0 if meanp > 200 else (np.abs(meanp-200))
        c[3] = 0 if max(self.linacfields) < 32 else (np.abs(max(self.linacfields)-32))
        c[4] = 0 if fhcfield < 35 else (np.abs(fhcfield-35))
        fitness = map(lambda x,y: (x*y)**2,c,w)
        return np.sqrt(np.sum(c))

def optfunc(args):
    fit = fitnessFunc(args)
    fitvalue = fit.calculateBeamParameters()
    fit.removeTempDirectory()
    return (fitvalue,)

# print optfunc([21,20,25,-25,25,-25,30,187,25,0])

def doOptimise():
    creator.create("FitnessMin", base.Fitness, weights=(1.0,))
    creator.create("Particle", list, fitness=creator.FitnessMin, speed=list,
    smin=None, smax=None, best=None)

    def generate(size, ranges, smin, smax):
        part = creator.Particle(random.uniform(a,b) for a,b in ranges)
        part.speed = [random.uniform(smin, smax) for _ in range(size)]
        part.smin = smin
        part.smax = smax
        return part

    def updateParticle(part, best, phi1, phi2):
        u1 = (random.uniform(0, phi1) for _ in range(len(part)))
        u2 = (random.uniform(0, phi2) for _ in range(len(part)))
        v_u1 = map(operator.mul, u1, map(operator.sub, part.best, part))
        v_u2 = map(operator.mul, u2, map(operator.sub, best, part))
        part.speed = list(map(operator.add, part.speed, map(operator.add, v_u1, v_u2)))
        for i, speed in enumerate(part.speed):
            if speed < part.smin:
                part.speed[i] = part.smin
            elif speed > part.smax:
                part.speed[i] = part.smax
        part[:] = list(map(operator.add, part, part.speed))

    toolbox = base.Toolbox()
    startranges = [[10, 32], [-30,30], [10, 32], [-30,30], [10, 32], [-30,30], [10, 32], [135,200], [10, 32], [-30,30]]
    toolbox.register("particle", generate, size=10, ranges=startranges, smin=-3, smax=3)
    toolbox.register("population", tools.initRepeat, list, toolbox.particle)
    toolbox.register("update", updateParticle, phi1=2.0, phi2=2.0)
    toolbox.register("evaluate", optfunc)

    pop = toolbox.population(n=5)
    stats = tools.Statistics(lambda ind: ind.fitness.values)
    stats.register("avg", np.mean)
    stats.register("std", np.std)
    stats.register("min", np.min)
    stats.register("max", np.max)

    logbook = tools.Logbook()
    logbook.header = ["gen", "evals"] + stats.fields

    GEN = 1000
    best = None

    import multiprocessing

    pool = multiprocessing.Pool()
    toolbox.register("map", pool.map)

    for g in range(GEN):
        for part in pop:
            part.fitness.values = toolbox.evaluate(part)
            if not part.best or part.best.fitness < part.fitness:
                part.best = creator.Particle(part)
                part.best.fitness.values = part.fitness.values
            if not best or best.fitness < part.fitness:
                best = creator.Particle(part)
                best.fitness.values = part.fitness.values
        for part in pop:
            toolbox.update(part, best)

        # Gather all the fitnesses in one list and print the stats
        logbook.record(gen=g, evals=len(pop), **stats.compile(pop))
        print(logbook.stream)

    return pop, logbook, best

if __name__ == '__main__':
    doOptimise()
