import numpy as np
import os, sys
sys.path.append(os.path.abspath(__file__+'/../../../../'))
from SimulationFramework.Framework import *
import FitnessFunc_Longitudinal_elegant as FF
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
opt = optimiser()
import shutil
import uuid

''' Run the injector part once if only optimising post-injector parameters'''
def create_base_files(scaling):
    framework = Framework('basefiles_'+str(scaling), overwrite=True)
    framework.loadSettings('Lattices/clara400_v12_v3.def')
    if not os.name == 'nt':
        framework.defineASTRACommand(['mpiexec','-np',str(3*scaling),'/opt/ASTRA/astra_MPICH2.sh'])
        framework.defineGeneratorCommand(['/opt/ASTRA/generator.sh'])
        framework.defineCSRTrackCommand(['/opt/OpenMPI-1.4.3/bin/mpiexec','-n',str(3*scaling),'/opt/CSRTrack/csrtrack_openmpi.sh'])
    framework.generator.number_of_particles = 2**(3*scaling)
    framework.track(files=['generator','injector400','S02'])

# for i in [3,4,5,6]:
#     create_base_files(i)
# exit()

def optfunc(args, dir=None, **kwargs):
    if dir == None:
        with FF.TemporaryDirectory(dir=os.getcwd()) as tmpdir:
            fit = FF.fitnessFunc(args, tmpdir, **kwargs)
            fitvalue = fit.calculateBeamParameters()
    else:
            fit = FF.fitnessFunc(args, dir, **kwargs)
            fitvalue = fit.calculateBeamParameters()
    return (fitvalue,)

framework = Framework('longitudinal_best', overwrite=False)
framework.loadSettings('Lattices/clara400_v12_v3_elegant.def')
injparameters = []
parameters = []
''' if including injector'''
injparameters.append(framework.getElement('CLA-HRG1-GUN-CAV', 'phase'))
injparameters.append(framework.getElement('CLA-HRG1-GUN-SOL', 'field_amplitude'))
injparameters.append(framework.getElement('CLA-L01-CAV', 'field_amplitude'))
injparameters.append(framework.getElement('CLA-L01-CAV', 'phase'))
injparameters.append(framework.getElement('CLA-L01-CAV-SOL-01', 'field_amplitude'))
injparameters.append(framework.getElement('CLA-L01-CAV-SOL-02', 'field_amplitude'))
''' always '''
parameters.append(framework.getElement('CLA-L02-CAV', 'field_amplitude'))
parameters.append(framework.getElement('CLA-L02-CAV', 'phase'))
parameters.append(framework.getElement('CLA-L03-CAV', 'field_amplitude'))
parameters.append(framework.getElement('CLA-L03-CAV', 'phase'))
parameters.append(framework.getElement('CLA-L4H-CAV', 'field_amplitude'))
parameters.append(framework.getElement('CLA-L4H-CAV', 'phase'))
parameters.append(framework.getElement('CLA-L04-CAV', 'field_amplitude'))
parameters.append(framework.getElement('CLA-L04-CAV', 'phase'))
parameters.append(framework.fileSettings['POSTINJ']['groups']['bunch_compressor']['dipoleangle'])
best = injparameters + parameters

#######################
post_inj = True
#######################

results = []
# with open('CLARA_longitudinal_best_solutions_simplex_elegant.csv.tmp', 'r') as csvfile:
#   reader = csv.reader(csvfile, delimiter=',', quoting=csv.QUOTE_NONNUMERIC)
#   for row in reader:
#     results.append(row)
# best = results[0][-len(parameters):]
best = [2.6338457327938296e7,-23.9868215448966,2.581910905052696e7,-7.618916138788988,2.43070395756709e7,188.3521131983386,2.7944819565259825e7,43.7590747875747,-0.1278008605127734]
if not post_inj:
    best = injparameters + best
print 'starting values = ', best

# fit = FF.fitnessFunc(best, os.getcwd()+'/test_3', scaling=3, overwrite=True, verbose=True, summary=False)
# print fit.calculateBeamParameters()
# exit()


# startranges = [[10, 32], [-40,40], [10, 32], [-40,40], [10, 50], [135,200],
#                [70, 100], [-40,40], [70, 100], [-40,40], [70, 100], [-40,40], [0.8,0.15]
#               ]
startranges = [[0.95*i, 1.05*i] if abs(i) > 0 else [-20,20] for i in best]
print 'startranges = ', startranges
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
    toolbox.register("evaluate", optfunc, scaling=5, post_injector=post_inj)
else:
    toolbox.register("evaluate", optfunc, scaling=5, post_injector=post_inj)

toolbox.register("mate", tools.cxBlend, alpha=0.2)
toolbox.register("mutate", tools.mutGaussian, mu=0, sigma=3, indpb=0.3)
toolbox.register("select", tools.selTournament, tournsize=3)


if __name__ == "__main__":
    global hof
    random.seed(64)

    if not post_inj:
        out = open('best_solutions_longitudinal_GA_elegant_Inj.csv','w', buffering=0)
        hoffile='CLARA_HOF_longitudinal_Elegant_Inj.csv'
    else:
        out = open('best_solutions_longitudinal_GA_elegant.csv','w', buffering=0)
        hoffile='CLARA_HOF_longitudinal_Elegant.csv'
    FF.csv_out = csv.writer(out)

    # Process Pool of 4 workers
    if not os.name == 'nt':
        pool = multiprocessing.Pool(processes=8)
    else:
        pool = multiprocessing.Pool(processes=3)
    toolbox.register("map", pool.map)

    if not os.name == 'nt':
        pop = toolbox.population(n=24)
    else:
        pop = toolbox.population(n=6)
    hof = tools.HallOfFame(10)
    stats = tools.Statistics(lambda ind: ind.fitness.values)
    stats.register("avg", np.mean)
    stats.register("std", np.std)
    stats.register("min", np.min)
    stats.register("max", np.max)

    pop, logbook = opt.eaSimple(pop, toolbox, cxpb=0.5, mutpb=0.2, ngen=50,
                            stats=stats, halloffame=hof, verbose=True,
                            hoffile=hoffile)

    pool.close()
    # print 'pop = ', pop
    print logbook
    print hof

    # try:
    print 'best fitness = ', optfunc(hof[0], dir=os.getcwd()+'/CLARA_best_longitudinal_elegant', scaling=6, overwrite=True, verbose=True, summary=True, post_injector=post_inj)
    #     with open('CLARA_best_longitudinal_elegant/CLARA_longitudinal_best_solutions_elegant.csv','wb') as out:
    #         csv_out=csv.writer(out)
    #         for row in hof:
    #             csv_out.writerow(row)
    #     with open('CLARA_best_longitudinal_elegant/CLARA_longitudinal_best_stats_elegant.csv','wb') as out:
    #         csv_out=csv.writer(out)
    #         for row in stats:
    #             csv_out.writerow(row)
    # except:
    #     with open('CLARA_longitudinal_best_solutions_elegant.csv.tmp','wb') as out:
    #         csv_out=csv.writer(out)
    #         for row in hof:
    #             csv_out.writerow(row)
