import os, sys
# sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname( os.path.abspath(__file__)))))
sys.path.append(os.path.abspath(__file__+'/../../../../'))
sys.path.append(os.path.abspath(__file__+'/../../../../../'))
from SimulationFramework.Framework import *
from SimulationFramework.Modules.constraints import *
import numpy as np
import SimulationFramework.Modules.read_twiss_file as rtf
import SimulationFramework.Modules.read_beam_file as rbf
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

def create_base_files(scaling):
    framework = Framework('long_240_'+str(scaling), overwrite=True)
    framework.loadSettings('Lattices/clara400_v12_v3.def')
    if not os.name == 'nt':
        framework.defineASTRACommand(['mpiexec','-np',str(3*scaling),'/opt/ASTRA/astra_MPICH2.sh'])
        framework.defineGeneratorCommand(['/opt/ASTRA/generator.sh'])
    framework.generator.number_of_particles = 2 ** (3 * scaling)
    results = []
    with open('CLARA_best_longitudinal_elegant\CLARA_longitudinal_best_solutions_elegant.csv', 'r') as csvfile:
        reader = csv.reader(csvfile, delimiter=',', quoting=csv.QUOTE_NONNUMERIC)
        for row in reader:
            results.append(row)
    best = results[0][0:9]
    framework.modifyElement('CLA-L02-CAV', 'field_amplitude', abs(best[0]))
    framework.modifyElement('CLA-L02-CAV', 'phase', best[1])
    framework.modifyElement('CLA-L03-CAV', 'field_amplitude', abs(best[2]))
    framework.modifyElement('CLA-L03-CAV', 'phase', best[3])
    framework.modifyElement('CLA-L4H-CAV', 'field_amplitude', abs(best[4]))
    framework.modifyElement('CLA-L4H-CAV', 'phase', best[5])
    framework.modifyElement('CLA-L04-CAV', 'field_amplitude', abs(best[6]))
    framework.modifyElement('CLA-L04-CAV', 'phase', best[7])
    framework['bunch_compressor'].set_angle(abs(best[8]))
    print best
    allbest = []
    with open('transverse_best_long_240_solutions.csv.tmp', 'r') as infile:
        reader = csv.reader(infile, quoting=csv.QUOTE_NONE, skipinitialspace=True)
        for row in reader:
            allbest.append(row)
    best = map(lambda x: float(x), allbest[0])
    parameters = list(best)
    framework.setElementType('quadrupole', 'k1', parameters)
    # framework.track(files=['generator','injector400','S02','L02','S03','L03','S04','L4H','S05','VBC','S06','L04','S07'])
    framework.track(files=['L02', 'S03', 'L03', 'S04', 'L4H', 'S05', 'VBC', 'S06', 'L04', 'S07'])

# for i in [4,5,6]:
#    create_base_files(i)
# exit()

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

framework = Framework('twiss_temp', overwrite=False)
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
parameters.append(framework.fileSettings['VBC']['groups']['bunch_compressor']['dipoleangle'])
best = injparameters + parameters

results = []
with open('CLARA_best_longitudinal_elegant\CLARA_longitudinal_best_solutions_elegant.csv', 'r') as csvfile:
  reader = csv.reader(csvfile, delimiter=',', quoting=csv.QUOTE_NONNUMERIC)
  for row in reader:
    results.append(row)
best = results[0][0:len(parameters)]
framework.modifyElement('CLA-L02-CAV', 'field_amplitude', abs(best[0]))
framework.modifyElement('CLA-L02-CAV', 'phase', best[1])
framework.modifyElement('CLA-L03-CAV', 'field_amplitude', abs(best[2]))
framework.modifyElement('CLA-L03-CAV', 'phase', best[3])
framework.modifyElement('CLA-L4H-CAV', 'field_amplitude', abs(best[4]))
framework.modifyElement('CLA-L4H-CAV', 'phase', best[5])
framework.modifyElement('CLA-L04-CAV', 'field_amplitude', abs(best[6]))
framework.modifyElement('CLA-L04-CAV', 'phase', best[7])
framework['bunch_compressor'].set_angle(abs(best[8]))
parameters = framework.getElementType('quadrupole','k1')
print 'parameters = ', parameters
best = parameters
allbest = []
with open('transverse_best_long_240_solutions.csv.tmp', 'r') as infile:
    reader = csv.reader(infile, quoting=csv.QUOTE_NONE, skipinitialspace=True)
    for row in reader:
        allbest.append(row)
best = map(lambda x: float(x), allbest[0])
parameters = list(best)
framework.setElementType('quadrupole', 'k1', parameters)
scaling=4
class fitnessFunc():

    def __init__(self, args, tempdir, scaling=4, overwrite=True, verbose=False, summary=False):
        self.cons = constraintsClass()
        self.beam = rbf.beam()
        self.twiss = rtf.twiss()
        self.tmpdir = tempdir
        self.verbose = verbose
        self.summary = summary
        self.parameters = list(args)
        self.dirname = os.path.basename(self.tmpdir)
        self.framework = Framework(self.dirname, clean=False)
        self.scaling=5
        self.framework.loadSettings('Lattices/clara400_v12_v3_elegant_jkj.def')
        with open('CLARA_best_longitudinal_elegant\CLARA_longitudinal_best_solutions_elegant.csv', 'r') as csvfile:
            reader = csv.reader(csvfile, delimiter=',', quoting=csv.QUOTE_NONNUMERIC)
            for row in reader:
                results.append(row)
        best = results[0][0:len(parameters)]
        allbest = []
        with open('transverse_best_long_240_solutions.csv.tmp', 'r') as infile:
            reader = csv.reader(infile, quoting=csv.QUOTE_NONE, skipinitialspace=True)
            for row in reader:
                allbest.append(row)
        quadbest = map(lambda x: float(x), allbest[0])
        quadparameters = list(quadbest)
        self.framework.setElementType('quadrupole', 'k1', quadparameters)
        self.framework.modifyElement('CLA-L02-CAV', 'field_amplitude', abs(best[0]))
        self.framework.modifyElement('CLA-L02-CAV', 'phase', best[1])
        self.framework.modifyElement('CLA-L03-CAV', 'field_amplitude', abs(best[2]))
        self.framework.modifyElement('CLA-L03-CAV', 'phase', best[3])
        self.framework.modifyElement('CLA-L4H-CAV', 'field_amplitude', abs(best[4]))
        self.framework.modifyElement('CLA-L4H-CAV', 'phase', best[5])
        self.framework.modifyElement('CLA-L04-CAV', 'field_amplitude', abs(best[6]))
        self.framework.modifyElement('CLA-L04-CAV', 'phase', best[7])
        self.framework['bunch_compressor'].set_angle(abs(best[8]))
        if not os.name == 'nt':
            self.framework.defineASTRACommand(['mpiexec','-np',str(3*scaling),'/opt/ASTRA/astra_MPICH2.sh'])
            self.framework.defineGeneratorCommand(['/opt/ASTRA/generator.sh'])
            self.framework.defineCSRTrackCommand(['/opt/OpenMPI-1.4.3/bin/mpiexec','-n',str(3*scaling),'/opt/CSRTrack/csrtrack_openmpi.sh'])
            # self.framework.generator.number_of_particles = 2**(3*scaling)
        else:
            self.framework.generator.number_of_particles = 2**(3*3)
        self.framework.defineElegantCommand(['elegant'])
        self.framework.setElementType('quadrupole','k1', self.parameters)

    def between(self, value, minvalue, maxvalue, absolute=True):
        if absolute:
            result = max([minvalue,min([maxvalue,abs(value)])])
        else:
            result = np.sign(value)*max([minvalue,min([maxvalue,abs(value)])])
        return result

    def calculateBeamParameters(self):
        twiss = self.twiss
        self.post_injector = True
        # try:
        if self.post_injector:
            startS = self.framework['POSTINJ'].startObject['position_start'][2]
            self.framework['POSTINJ'].file_block['input']['prefix'] = '../basefiles_'+str(self.scaling)+'/'
            # self.framework.modifyElement('CLA-L02-CAV', 'field_amplitude', abs(best[0]))
            # self.framework.modifyElement('CLA-L02-CAV', 'phase', best[1])
            # self.framework.modifyElement('CLA-L03-CAV', 'field_amplitude', abs(best[2]))
            # self.framework.modifyElement('CLA-L03-CAV', 'phase', best[3])
            # self.framework.modifyElement('CLA-L4H-CAV', 'field_amplitude', abs(best[4]))
            # self.framework.modifyElement('CLA-L4H-CAV', 'phase', best[5])
            # self.framework.modifyElement('CLA-L04-CAV', 'field_amplitude', abs(best[6]))
            # self.framework.modifyElement('CLA-L04-CAV', 'phase', best[7])
            # self.framework['bunch_compressor'].set_angle(abs(best[8]))
            self.framework.track(startfile='POSTINJ')
        else:
            self.framework.track()
        constraintsList = {}
        constraintsListQuads = {
            'max_k': {'type': 'lessthan', 'value': [abs(p) for p in self.parameters], 'limit': 2.5, 'weight': 10},

        }
        twiss.reset_dicts()
        twiss.read_sdds_file(self.dirname + '/' + 'POSTINJ.twi')
        twiss.read_sdds_file(self.dirname + '/' + 'POSTINJ.sig')
        twiss['s'] += startS
        # constraintsList = merge_two_dicts(constraintsList, constraintsListQuads)
        # twiss.read_astra_emit_files( [ self.dirname+'/'+n+'.Zemit.001' for n in self.framework.fileSettings.keys() if self.framework.fileSettings[n]['code'].upper() == 'ASTRA'] )
        constraintsListSigmas = {
            'max_xrms': {'type': 'lessthan', 'value': 1e3*twiss['Sx'], 'limit': 1, 'weight': 10},
            'max_yrms': {'type': 'lessthan', 'value': 1e3*twiss['Sy'], 'limit': 1, 'weight': 10},
            # 'max_betax': {'type': 'lessthan', 'value': twiss['betax'], 'limit': 100, 'weight': 100},
            'min_xrms': {'type': 'greaterthan', 'value': 1e3*twiss['Sx'], 'limit': 0.1, 'weight': 10},
            'min_yrms': {'type': 'greaterthan', 'value': 1e3*twiss['Sy'], 'limit': 0.1, 'weight': 10},
            'last_exn': {'type': 'lessthan', 'value': 1e6*twiss['enx'][-1], 'limit': 0.6, 'weight': 45},
            'last_eyn': {'type': 'lessthan', 'value': 1e6*twiss['eny'][-1], 'limit': 0.6, 'weight': 1},
        }
        constraintsList = merge_two_dicts(constraintsList, constraintsListSigmas)
        # twiss.read_astra_emit_files(self.dirname+'/S07.Zemit.001')
        tdc_position = self.framework['CLA-S07-TDC-01-R']['position_start'][2]
        tdc_screen_position = self.framework['CLA-S07-DIA-SCR-03-W']['position_start'][2]
        dechirper_position = self.framework['CLA-S07-DCP-01']['position_start'][2]
        constraintsListS07 = {
            'tdc_phase_advance': {'type': 'equalto', 'value': twiss.interpolate(tdc_screen_position,'psiy', index='s') - twiss.interpolate(tdc_position,'psiy', index='s'), 'limit': 0.25, 'weight': 1},
            'tdc_screen_beta_y': {'type': 'greaterthan', 'value': twiss.extract_values('betay', tdc_position, tdc_screen_position), 'limit': 5, 'weight': 1},
            'tdc_screen_beta_x': {'type': 'lessthan', 'value': twiss.extract_values('betax', tdc_position, tdc_screen_position), 'limit': 20, 'weight': 100},
            'dechirper_sigma_x': {'type': 'equalto', 'value': 1e3*twiss.interpolate(dechirper_position, 'Sx', index='s'), 'limit': 0.1, 'weight': 10},
            'dechirper_sigma_y': {'type': 'equalto', 'value': 1e3*twiss.interpolate(dechirper_position, 'Sy', index='s'), 'limit': 0.1, 'weight': 10},
            'dechirper_sigma_xy': {'type': 'equalto', 'value': 1e3*twiss.interpolate(dechirper_position, 'Sy', index='s'), 'limit': 1e3*twiss.interpolate(dechirper_position, 'Sx', index='s'), 'weight': 20},
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

allbest=[]
with open('transverse_best_long_240_solutions.csv.tmp','r') as infile:
    reader = csv.reader(infile, quoting=csv.QUOTE_NONE, skipinitialspace=True)
    for row in reader:
        allbest.append(row)
best = map(lambda x: float(x), allbest[0])
# best = parameters
# print 'starting values = ', best
# print optfunc(best, dir=os.getcwd()+'/test_transverse', scaling=6, overwrite=True, verbose=True, summary=False)
# exit()


# startranges = [[10, 32], [-40,40], [10, 32], [-40,40], [10, 32], [-40,40], [10, 32], [135,200], [10, 32], [-40,40], [0.8,0.15]]
startranges = [[0.9*i, 1.1*i] if abs(i) > 0 else [-0.1,0.1] for i in best]
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
    toolbox.register("evaluate", optfunc, scaling=3)
else:
    toolbox.register("evaluate", optfunc, scaling=4)
toolbox.register("mate", tools.cxBlend, alpha=0.2)
toolbox.register("mutate", tools.mutGaussian, mu=0, sigma=3, indpb=0.3)
toolbox.register("select", tools.selTournament, tournsize=3)


if __name__ == "__main__":
    random.seed(64)

    # Process Pool of 4 workers
    if not os.name == 'nt':
        pool = multiprocessing.Pool(processes=6)
    else:
        pool = multiprocessing.Pool(processes=3)
    toolbox.register("map", pool.map)
    # toolbox.register("map", futures.map)

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

    pop, logbook = algorithms.eaSimple(pop, toolbox, cxpb=0.5, mutpb=0.3, ngen=30,
                            stats=stats, halloffame=hof)

    # print 'pop = ', pop
    print logbook
    print hof

    try:
        print 'best fitness = ', optfunc(hof[0][0:31], dir=os.getcwd()+'/transverse_best_long_240', scaling=6, overwrite=True, verbose=True, summary=True)
        with open('transverse_best_long_240/transverse_best_solutions.csv','wb') as out:
            csv_out=csv.writer(out)
            for row in hof:
                csv_out.writerow(row)
    except:
        with open('transverse_best_long_240_solutions.csv.tmp','wb') as out:
            csv_out=csv.writer(out)
            for row in hof:
                csv_out.writerow(row)
    pool.close()
