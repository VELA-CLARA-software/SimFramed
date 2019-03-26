import os, sys
sys.path.append('../../../../')
import SimulationFramework.Framework as fw
import numpy as np
from SimulationFramework.Modules.constraints import *
import SimulationFramework.Modules.read_twiss_file as rtf
import SimulationFramework.Modules.read_beam_file as rbf
from SimulationFramework.Modules.nelder_mead import nelder_mead
from functools import partial
from scipy.optimize import minimize
import csv
from ruamel import yaml
import shutil
import uuid

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

framework = fw.Framework(None)
framework.loadSettings('Lattices/clara400_v12_v3.def')
parameters = framework.getElementType('quadrupole','k1l')
names = framework.getElementType('quadrupole','objectname')
# print 'parameters = ', parameters
best = parameters
scaling=4

best = [-0.22789817,  0.04384427,  0.07237042,  0.11319594,  0.11633546, -0.10393746,
  0.09247306,  0.03135896, -0.06080841, -0.04500804,  0.02695322,  0.0206167,
  0.03058594, -0.04103264, -0.12178037,  0.30441347,  0.22639786, -0.08582796,
  0.37171019, -0.13076231, -0.33075536, -0.10981188,  0.43603262,  0.01990002,
  0.27849027,  0.37018414, -0.1794828,   0.03491095, -0.17280846]

with open('best_changes.yaml', 'r') as infile:
    data = dict(yaml.load(infile, Loader=yaml.UnsafeLoader))
    best = [data[n]['k1l'] for n in names[:-2]]

# print best
# exit()

class fitnessFunc():

    def __init__(self, args, tempdir, scaling=4, overwrite=True, verbose=False, summary=False):
        self.cons = constraintsClass()
        self.beam = rbf.beam()
        self.twiss = rtf.twiss()
        self.tmpdir = tempdir
        self.verbose = verbose
        self.summary = summary
        self.parameters = list(args)
        self.dirname = os.path.relpath(self.tmpdir)
        self.framework = fw.Framework(self.dirname, clean=True, verbose=False)
        self.framework.loadSettings('Lattices/clara400_v12_v3_elegant_jkj.def')
        self.framework.change_Lattice_Code('All','elegant')
        if not os.name == 'nt':
            self.framework.defineASTRACommand(['mpiexec','-np',str(3*scaling),'/opt/ASTRA/astra_MPICH2.sh'])
            self.framework.defineGeneratorCommand(['/opt/ASTRA/generator.sh'])
            self.framework.defineCSRTrackCommand(['/opt/OpenMPI-1.4.3/bin/mpiexec','-n',str(3*scaling),'/opt/CSRTrack/csrtrack_openmpi.sh'])
            self.framework.generator.number_of_particles = 2**(3*scaling)
        else:
            self.framework.generator.number_of_particles = 2**(3*3)
        self.framework.defineElegantCommand(['elegant'])
        self.framework.setElementType('quadrupole','k1l', self.parameters + [0,0])

    def between(self, value, minvalue, maxvalue, absolute=True):
        if absolute:
            result = max([minvalue,min([maxvalue,abs(value)])])
        else:
            result = np.sign(value)*max([minvalue,min([maxvalue,abs(value)])])
        return result

    def calculateBeamParameters(self):
        # try:
            twiss = self.twiss
            self.framework['POSTINJ'].prefix = '../../../basefiles_4/'
            # print 'before tracking '
            self.framework.track(startfile='POSTINJ')
            # print 'after tracking '
            constraintsList = {}
            quadkls = self.framework.getElementType('quadrupole','k1l')
            quadlengths = self.framework.getElementType('quadrupole','length')
            quadnames = self.framework.getElementType('quadrupole','objectname')
            constraintsListQuads = {
                'max_k': {'type': 'lessthan', 'value': [abs(k/l) for k, l in zip(quadkls, quadlengths)], 'limit': 2.5, 'weight': 10},

            }
            constraintsList = merge_two_dicts(constraintsList, constraintsListQuads)

            twiss.read_elegant_twiss_files( [ self.dirname+'/'+n+'.twi' for n in ['POSTINJ']])
            constraintsListSigmas = {
                'max_xrms': {'type': 'lessthan', 'value': 1e3*twiss['sigma_x'], 'limit': 1, 'weight': 10},
                'max_yrms': {'type': 'lessthan', 'value': 1e3*twiss['sigma_y'], 'limit': 1, 'weight': 10},
                'min_xrms': {'type': 'greaterthan', 'value': 1e3*twiss['sigma_x'], 'limit': 0.1, 'weight': 10},
                'min_yrms': {'type': 'greaterthan', 'value': 1e3*twiss['sigma_y'], 'limit': 0.1, 'weight': 10},
                'last_exn': {'type': 'lessthan', 'value': 1e6*twiss['enx'][-1], 'limit': 0.6, 'weight': 1},
                'last_eyn': {'type': 'lessthan', 'value': 1e6*twiss['eny'][-1], 'limit': 0.6, 'weight': 1},
            }
            constraintsList = merge_two_dicts(constraintsList, constraintsListSigmas)

            twiss.read_elegant_twiss_files(self.dirname+'/POSTINJ.twi')
            tdc_position = self.framework['CLA-S07-TDC-01-R']['position_start'][2]
            tdc_screen_position = self.framework['CLA-S07-DIA-SCR-03-W']['position_start'][2]
            dechirper_position = self.framework['CLA-S07-DCP-01']['position_start'][2]
            constraintsListS07 = {
                'tdc_phase_advance': {'type': 'equalto', 'value': twiss.interpolate(tdc_screen_position,'muy') - twiss.interpolate(tdc_position,'muy'), 'limit': 0.25, 'weight': 1},
                'tdc_screen_beta_y': {'type': 'greaterthan', 'value': twiss.extract_values('beta_y', tdc_position, tdc_screen_position), 'limit': 5, 'weight': 1},
                'dechirper_sigma_x': {'type': 'equalto', 'value': 1e3*twiss.interpolate(dechirper_position, 'sigma_x'), 'limit': 0.1, 'weight': 10},
                'dechirper_sigma_y': {'type': 'equalto', 'value': 1e3*twiss.interpolate(dechirper_position, 'sigma_y'), 'limit': 0.1, 'weight': 10},
                'dechirper_sigma_xy': {'type': 'equalto', 'value': 1e3*twiss.interpolate(dechirper_position, 'sigma_y') - 1e3*twiss.interpolate(dechirper_position, 'sigma_x'), 'limit': 0.0, 'weight': 20},
            }
            constraintsList = merge_two_dicts(constraintsList, constraintsListS07)
            fitness = self.cons.constraints(constraintsList)
            if self.verbose:
                print self.cons.constraintsList(constraintsList)
            if self.summary:
                self.astra.createHDF5Summary(reference='Transverse_GA')
            return fitness
        # except:
        #     return 1e6

def optfunc(args, dir=None, **kwargs):
    global simplex_iteration, bestfit
    dir='simplex_transverse/simplex_iteration_'+str(simplex_iteration)
    if dir == None:
        with TemporaryDirectory(dir=os.getcwd()) as tmpdir:
            fit = fitnessFunc(args, tmpdir, **kwargs)
            fitvalue = fit.calculateBeamParameters()
    else:
            fit = fitnessFunc(args, dir, **kwargs)
            fitvalue = fit.calculateBeamParameters()
    fit.framework.save_changes_file(filename=fit.framework.subdirectory+'/changes.yaml', type='quadrupole')
    print 'fitvalue[', simplex_iteration, '] = ', fitvalue
    saveState(args, simplex_iteration, fitvalue)
    simplex_iteration += 1
    if fitvalue < bestfit:
        print '!!!!!!  New best = ', fitvalue
        fit.framework.save_changes_file(filename='best_changes.yaml', type='quadrupole')
        bestfit = fitvalue
    return fitvalue

savefile = open('simplex_transverse/best_solutions_running.csv','a')

def saveState(args, n, fitness):
    global savefile
    csv_out=csv.writer(savefile)
    args=list(args)
    args.append(n)
    args.append(fitness)
    csv_out.writerow(args)

simplexfunc = partial(optfunc, verbose=False)
bestfit = 1e26
with open('simplex_transverse/best_solutions_running.csv','w') as out:
    simplex_iteration = 0

# res = minimize(simplexfunc, best, method='nelder-mead', options={'disp': True, 'maxiter': 300, 'adaptive': True})
# print list(res.x)

res = nelder_mead(optfunc, np.array(best), step=0.05, max_iter=300, no_improv_break=100)
print res
