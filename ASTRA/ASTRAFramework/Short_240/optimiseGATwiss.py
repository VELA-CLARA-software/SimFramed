from ASTRAInjector import *
from CSRTrack import *
import numpy as np
from constraints import *
import os
import tempfile
import copy
import read_twiss_file as rtf
import read_beam_file as raf
from deap import algorithms
from deap import base
from deap import creator
from deap import tools
import multiprocessing
from scoop import futures
import operator
import random
from shutil import copyfile
import csv
twiss = rtf.twiss()
beam = raf.beam()


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

    def __init__(self, args, tempdir, npart=4096, ncpu=4, overwrite=True, verbose=False, summary=False):
        self.tmpdir = tempdir
        self.verbose = verbose
        self.parameters = list(args)
        self.scgrid = int(round(npart**(1./3.)))
        # print [abs(p) for p in self.parameters]
        self.dirname = os.path.basename(self.tmpdir)
        astra = ASTRAInjector(self.dirname, overwrite=overwrite)
        csrtrack = CSRTrack(self.dirname, overwrite=overwrite)
        if not os.name == 'nt':
            astra.defineASTRACommand(['mpiexec','-np',str(ncpu),'/opt/ASTRA/astra_MPICH2.sh'])
            csrtrack.defineCSRTrackCommand(['/opt/OpenMPI-1.4.3/bin/mpiexec','-n',str(ncpu),'/opt/CSRTrack/csrtrack_openmpi.sh'])
        else:
            astra.defineASTRACommand(['astra'])
            csrtrack.defineCSRTrackCommand(['CSRtrack_1.201.wic.exe'])
        astra.loadSettings('short_240_12b3.settings')
        scgrid = getGrids(npart)
        astra.globalSettings['SC_2D_Nrad'] = max([scgrid.gridSizes, 4])
        astra.globalSettings['SC_2D_Nlong'] = max([scgrid.gridSizes, 4])
        for scvar in ['SC_3D_Nxf', 'SC_3D_Nyf', 'SC_3D_Nzf']:
            astra.globalSettings[scvar] = scgrid.gridSizes
        astra.fileSettings['short_240.2']['quad_K'] = self.parameters[0:6]
        astra.fileSettings['short_240.3']['quad_K'] = self.parameters[6:14]
        astra.fileSettings['short_240.4']['quad_K'] = self.parameters[14:16]
        astra.fileSettings['short_240.5']['quad_K'] = self.parameters[16:]
        bcangle = float(astra.fileSettings['vb']['variable_bunch_compressor']['angle'])
        # print 'Creating Initial Distribution in folder:', self.tmpdir
        astra.createInitialDistribution(npart=npart, charge=250)
        # print 'Apply Settings in folder:', self.tmpdir
        astra.fileSettings['short_240.5']['starting_distribution'] = 'end.fmt2.astra'
        ''' Create ASTRA files based on settings '''
        astra.applySettings()
        ''' Run ASTRA upto VBC '''
        astra.runASTRAFiles(files=['short_240.1','short_240.2','short_240.3','short_240.4'])
        ''' Write Out the CSRTrack file based on the BC angle (assumed to be 0.105) '''
        csrtrack.writeCSRTrackFile('csrtrk.in', angle=bcangle, forces='projected')
        ''' Run CSRTrack'''
        csrtrack.runCSRTrackFile('csrtrk.in')
        ''' Convert CSRTrack output file back in to ASTRA format '''
        beam.convert_csrtrackfile_to_astrafile(self.dirname+'/'+'end.fmt2', self.dirname+'/'+'end.fmt2.astra')
        ''' Run the next section of the lattice in ASTRA, using the CSRTrack output as input '''
        astra.runASTRAFiles(files=['short_240.5'])
        self.cons = constraintsClass()
        if summary:
            astra.createHDF5Summary()

    def calculateTwissParameters(self):
        constraintsList = {}
        # LINAC 2 and 3
        constraintsListQuads = {
            'max_k': {'type': 'lessthan', 'value': [abs(p) for p in self.parameters], 'limit': 0.1, 'weight': 2},

        }
        constraintsList = merge_two_dicts(constraintsList, constraintsListQuads)
        twiss.read_astra_emit_files(self.dirname+'/short_240.2.Zemit.001')
        constraintsList2 = {
            'max_xrms_2': {'type': 'lessthan', 'value': 1e3*twiss['sigma_x'], 'limit': 1, 'weight': 10},
            'max_yrms_2': {'type': 'lessthan', 'value': 1e3*twiss['sigma_y'], 'limit': 1, 'weight': 10},
            'min_xrms_2': {'type': 'greaterthan', 'value': 1e3*twiss['sigma_x'], 'limit': 0.2, 'weight': 10},
            'min_yrms_2': {'type': 'greaterthan', 'value': 1e3*twiss['sigma_y'], 'limit': 0.2, 'weight': 10},
            'last_exn_2': {'type': 'lessthan', 'value': 1e6*twiss['enx'][-1], 'limit': 0.8, 'weight': 100},
            'last_eyn_2': {'type': 'lessthan', 'value': 1e6*twiss['eny'][-1], 'limit': 0.8, 'weight': 100},
            # 'beta_x_2': {'type': 'lessthan', 'value': twiss['beta_x'], 'limit': 50, 'weight': 10},
            # 'beta_y_2': {'type': 'lessthan', 'value': twiss['beta_y'], 'limit': 50, 'weight': 10},
        }
        constraintsList = merge_two_dicts(constraintsList, constraintsList2)
        # 4HC
        twiss.read_astra_emit_files(self.dirname+'/short_240.3.Zemit.001')
        constraintsList3 = {
            'max_xrms_3': {'type': 'lessthan', 'value': 1e3*twiss['sigma_x'], 'limit': 1, 'weight': 10},
            'max_yrms_3': {'type': 'lessthan', 'value': 1e3*twiss['sigma_y'], 'limit': 1, 'weight': 10},
            'min_xrms_3': {'type': 'greaterthan', 'value': 1e3*twiss['sigma_x'], 'limit': 0.2, 'weight': 20},
            'min_yrms_3': {'type': 'greaterthan', 'value': 1e3*twiss['sigma_y'], 'limit': 0.2, 'weight': 20},
            'last_exn_3': {'type': 'lessthan', 'value': 1e6*twiss['enx'][-1], 'limit': 0.8, 'weight': 100},
            'last_eyn_3': {'type': 'lessthan', 'value': 1e6*twiss['eny'][-1], 'limit': 0.8, 'weight': 100},
            # 'beta_x_3': {'type': 'lessthan', 'value': twiss['beta_x'], 'limit': 50, 'weight': 10},
            # 'beta_y_3': {'type': 'lessthan', 'value': twiss['beta_y'], 'limit': 50, 'weight': 10},
        }
        constraintsList = merge_two_dicts(constraintsList, constraintsList3)

        # VBC
        twiss.read_astra_emit_files(self.dirname+'/short_240.4.Zemit.001')
        constraintsList4 = {
            'min_xrms_4': {'type': 'greaterthan', 'value': 1e3*twiss['sigma_x'], 'limit': 0.2, 'weight': 50},
            'min_yrms_4': {'type': 'greaterthan', 'value': 1e3*twiss['sigma_y'], 'limit': 0.2, 'weight': 50},
            'last_yrms_4': {'type': 'greaterthan', 'value': 1e3*twiss['sigma_y'][-1], 'limit': 0.4, 'weight': 60},
            'last_exn_4': {'type': 'lessthan', 'value': 1e6*twiss['enx'][-1], 'limit': 0.8, 'weight': 100},
            'last_eyn_4': {'type': 'lessthan', 'value': 1e6*twiss['eny'][-1], 'limit': 0.8, 'weight': 100},
            # 'beta_x_4': {'type': 'lessthan', 'value': twiss['beta_x'], 'limit': 50, 'weight': 10},
            # 'beta_y_4': {'type': 'lessthan', 'value': twiss['beta_y'], 'limit': 50, 'weight': 10},
        }
        ''' This doesn't make much sense with CSRTRack being used, but still... '''
        # constraintsList = merge_two_dicts(constraintsList, constraintsList4)
        # LINAC 4 and TDC (40.8012m) with screen (46.4378m) and Dechirper (44.03m)
        twiss.read_astra_emit_files(self.dirname+'/short_240.5.Zemit.001')
        constraintsList5 = {
            'max_xrms_5': {'type': 'lessthan', 'value': 1e3*twiss['sigma_x'], 'limit': 1, 'weight': 10},
            'max_yrms_5': {'type': 'lessthan', 'value': 1e3*twiss['sigma_y'], 'limit': 1, 'weight': 10},
            'min_xrms_5': {'type': 'greaterthan', 'value': 1e3*twiss['sigma_x'], 'limit': 0.2, 'weight': 20},
            'min_yrms_5': {'type': 'greaterthan', 'value': 1e3*twiss['sigma_y'], 'limit': 0.2, 'weight': 20},
            # 'last_alpha_x_5': {'type': 'lessthan', 'value': abs(twiss['alpha_x'][-1]), 'limit': 2, 'weight': 10},
            # 'last_alpha_y_5': {'type': 'lessthan', 'value': abs(twiss['alpha_y'][-1]), 'limit': 2, 'weight': 10},
            # 'last_beta_x_5': {'type': 'lessthan', 'value': twiss['beta_x'], 'limit': 30, 'weight': 2.},
            # 'last_beta_y_5': {'type': 'lessthan', 'value': twiss['beta_y'], 'limit': 30, 'weight': 2.},
            'tdc_phase_advance': {'type': 'equalto', 'value': twiss.interpolate(46.4378,'muy') - twiss.interpolate(40.8012,'muy'), 'limit': 0.25, 'weight': 500},
            # 'tdc_beta_y_greaterthan': {'type': 'greaterthan', 'value': twiss.interpolate(40.8012, 'beta_y'), 'limit': 50, 'weight': 25},
            # 'tdc_beta_y_lassthan': {'type': 'lessthan', 'value': twiss.interpolate(40.8012, 'beta_y'), 'limit': 100, 'weight': 25},
            'tdc_screen_beta_y': {'type': 'greaterthan', 'value': twiss.extract_values('beta_y', 40.8012, 46.4378), 'limit': 5, 'weight': 50},
            # 'screen_beta_x': {'type': 'equalto', 'value': twiss.interpolate(46.4378, 'beta_x'), 'limit': 5, 'weight': 75},
            # 'screen_beta_y': {'type': 'equalto', 'value': twiss.interpolate(46.4378, 'beta_y'), 'limit': 5, 'weight': 75},
            'dechirper_sigma_x': {'type': 'lessthan', 'value': twiss.interpolate(44.03, 'sigma_x'), 'limit': 0.1, 'weight': 100},
            'dechirper_sigma_y': {'type': 'lessthan', 'value': twiss.interpolate(44.03, 'sigma_y'), 'limit': 0.1, 'weight': 100},
            'last_exn_5': {'type': 'lessthan', 'value': 1e6*twiss['enx'], 'limit': 0.8, 'weight': 200},
            'last_eyn_5': {'type': 'lessthan', 'value': 1e6*twiss['eny'], 'limit': 0.8, 'weight': 200},
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

# results = []

# with open('twiss_best/twiss_best_solutions.csv', 'r') as csvfile:
#   reader = csv.reader(csvfile, delimiter=',', quoting=csv.QUOTE_NONNUMERIC)
#   for row in reader:
#     results.append(row)
# best = results[0]
astra = ASTRAInjector('twiss_best', overwrite=False)
astra.loadSettings('short_240_12b3.settings')
parameters = []
parameters.append(astra.fileSettings['short_240.2']['quad_K'])
parameters.append(astra.fileSettings['short_240.3']['quad_K'])
parameters.append(astra.fileSettings['short_240.4']['quad_K'])
parameters.append(astra.fileSettings['short_240.5']['quad_K'])
best = [item for sublist in parameters for item in sublist]

# print optfunc(best, dir=os.getcwd()+'/twiss_best', npart=50000, ncpu=20, overwrite=False, verbose=False, summary=True)
# exit()

# best = [0 for x in best]

startranges = [[0.8*i, 1.2*i] if abs(i) > 0 else [-0.1,0.1] for i in best]
print 'Start Ranges = ', startranges
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
    toolbox.register("evaluate", optfunc, npart=2**(3*3))
else:
    toolbox.register("evaluate", optfunc, npart=2**(4*3))
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

    pop, logbook = algorithms.eaSimple(pop, toolbox, cxpb=0.5, mutpb=0.2, ngen=10,
                            stats=stats, halloffame=hof)

    # print 'pop = ', pop
    print logbook
    print hof

    try:
        with open('twiss_best/twiss_best_solutions.csv','wb') as out:
            csv_out=csv.writer(out)
            for row in hof:
                csv_out.writerow(row)
    except:
        with open('twiss_best/twiss_best_solutions.csv.tmp','wb') as out:
            csv_out=csv.writer(out)
            for row in hof:
                csv_out.writerow(row)
    pool.close()

    print 'best fitness = ', optfunc(hof[0], dir=os.getcwd()+'/twiss_best', npart=4096, ncpu=40)
