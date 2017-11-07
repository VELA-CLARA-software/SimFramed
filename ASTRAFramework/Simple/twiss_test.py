from ASTRAInjector import *
import numpy as np
from constraints import *
import os
import copy
import read_twiss_file as rtf
from scipy.optimize import minimize
import operator
import random

import csv

twiss = rtf.twiss()

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

    def __init__(self, args, tempdir, npart=100, ncpu=10, overwrite=True, verbose=False):
        self.tmpdir = tempdir
        self.verbose = verbose
        self.parameters = list(args)
        self.dirname = os.path.basename(self.tmpdir)
        astra = ASTRAInjector(self.dirname, overwrite=overwrite)
        if not os.name == 'nt':
            astra.defineASTRACommand(['mpiexec','-np',str(ncpu),'/opt/ASTRA/astra_MPICH2_Quiet.sh'])
        astra.loadSettings('short_240_0955.settings')
        astra.fileSettings['test.in.127']['quad_K'] = self.parameters[0:2]
        astra.fileSettings['test.in.128']['quad_K'] = self.parameters[2:]
        astra.createInitialDistribution(npart=npart, charge=250)
        astra.applySettings()
        astra.runASTRAFiles()
        # ft = feltools(self.dirname)
        # sddsfile = ft.convertToSDDS('test.in.128.4929.128')
        self.cons = constraintsClass()

    def calculateTwissParameters(self):
        # twiss.read_astra_emit_files('test.in.126.Zemit.126')
        # constraintsList126 = {
        #     'max_xrms_126': {'type': 'lessthan', 'value': max(twiss['sigma_x']), 'limit': 1, 'weight': 10},
        #     'max_yrms_126': {'type': 'lessthan', 'value': max(twiss['sigma_y']), 'limit': 1, 'weight': 10},
        #     'min_xrms_126': {'type': 'greaterthan', 'value': min(twiss['sigma_x']), 'limit': 0.2, 'weight': 10},
        #     'min_yrms_126': {'type': 'greaterthan', 'value': min(twiss['sigma_y']), 'limit': 0.2, 'weight': 10},
        # }
        twiss.read_astra_emit_files(self.dirname+'/test.in.127.Zemit.127')
        constraintsList127 = {
            'min_xrms_127': {'type': 'greaterthan', 'value': 1e3*min(twiss['sigma_x']), 'limit': 0.1, 'weight': 30},
            'min_yrms_127': {'type': 'greaterthan', 'value': 1e3*min(twiss['sigma_y']), 'limit': 0.2, 'weight': 50},
            'last_yrms_127': {'type': 'greaterthan', 'value': 1e3*twiss['sigma_y'][-1], 'limit': 0.4, 'weight': 30},
            'last_xrms_127': {'type': 'lessthan', 'value': 1e3*twiss['sigma_x'][-1], 'limit': 0.2, 'weight': 30},
            'last_exn_127': {'type': 'lessthan', 'value': 1e6*twiss['enx'][-1], 'limit': 0.8, 'weight': 10},
            'last_eyn_127': {'type': 'lessthan', 'value': 1e6*twiss['eny'][-1], 'limit': 0.8, 'weight': 10},
        }
        twiss.read_astra_emit_files(self.dirname+'/test.in.128.Zemit.128')
        constraintsList128 = {
            'min_xrms_128': {'type': 'greaterthan', 'value': 1e3*min(twiss['sigma_x']), 'limit': 0.1, 'weight': 20},
            'min_yrms_128': {'type': 'greaterthan', 'value': 1e3*min(twiss['sigma_y']), 'limit': 0.2, 'weight': 20},
            'last_alpha_x_128': {'type': 'lessthan', 'value': abs(twiss['alpha_x'][-1]), 'limit': 5, 'weight': 10},
            'last_alpha_y_128': {'type': 'lessthan', 'value': abs(twiss['alpha_y'][-1]), 'limit': 5, 'weight': 10},
            'last_beta_x_128': {'type': 'lessthan', 'value': twiss['beta_x'][-1], 'limit': 50, 'weight': 50},
            'last_beta_y_128': {'type': 'lessthan', 'value': twiss['beta_y'][-1], 'limit': 50, 'weight': 50},
        }
        constraintsList = merge_two_dicts(constraintsList127, constraintsList128)
        fitness = self.cons.constraints(constraintsList)
        if self.verbose:
            print self.cons.constraintsList(constraintsList)
        return fitness

def optfunc(args, dir=None, **kwargs):
    if dir == None:
        with TemporaryDirectory(dir=os.getcwd()) as tmpdir:
            fit = fitnessFunc(args, tmpdir, **kwargs)
            fitvalue = fit.calculateTwissParameters()
    else:
            fit = fitnessFunc(args, dir, **kwargs)
            fitvalue = fit.calculateTwissParameters()
    return fitvalue

if not os.name == 'nt':
    os.chdir('/home/jkj62/ASTRAFramework/Simple')

best = [
    1.76571, -0.334896,
    -2.38819, 3.60297, 5., -4., -3.29565, 3.38181, 1.52921, 2.22559, -2.95347, -1.71192, -5.98555, 7.72644, -5.24599
]
# print optfunc(best, dir=os.getcwd()+'/match_1117', npart=1000, ncpu=20, overwrite=False, verbose=True)
# exit()

startranges = [[0.8*i, 1.4*i] for i in best]


def createSimplex():
    simplex = [copy.copy(best) for i in range(len(best)+1)]
    for i in range(len(best)):
        lower, upper = startranges[i]
        r = random.uniform(lower, upper)
        simplex[i][i] = r
    return simplex

if __name__ == "__main__":
    random.seed(64)

    startsimplex = createSimplex()
    # print startsimplex
    # exit()
    res = minimize(optfunc, best, method='Nelder-Mead', tol=1e-6, options={'initial_simplex': startsimplex, 'disp': True})
    print res
    print 'best fitness = ', optfunc(res.x, dir=os.getcwd()+'/match', npart=10000, ncpu=40)
