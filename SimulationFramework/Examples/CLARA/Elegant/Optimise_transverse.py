import os, sys
sys.path.append('../../../../')
import SimulationFramework.Framework as fw
import numpy as np
from SimulationFramework.Modules.nelder_mead import nelder_mead
import SimulationFramework.Examples.CLARA.Elegant.runElegant as runEle
from SimulationFramework.Modules.merge_two_dicts import merge_two_dicts
from SimulationFramework.Modules.constraints import *
from functools import partial
from scipy.optimize import minimize
import csv
from ruamel import yaml
import shutil
import uuid

class Optimise_transverse(runEle.fitnessFunc):

    def __init__(self, lattice='Lattices/clara400_v12_v3.def', scaling=5):
        super(Optimise_transverse, self).__init__()
        self.lattice_file = lattice
        self.framework = fw.Framework(None)
        self.framework.loadSettings(self.lattice_file)
        self.parameter_names = [q for q in self.framework.getElementType('quadrupole','objectname')]
        self.parameters = [[q, 'k1l'] for q in self.framework.getElementType('quadrupole','objectname')]
        self.cons = constraintsClass()
        self.changes = None
        self.resultsDict = {}
        # ******************************************************************************
        CLARA_dir = os.path.relpath(__file__+'/../../')
        self.POST_INJECTOR = True
        CREATE_BASE_FILES = False
        if self.POST_INJECTOR and CREATE_BASE_FILES:
            self.optfunc = partial(self.OptimisingFunction, scaling=scaling, post_injector=self.POST_INJECTOR, basefiles='../CLARA/basefiles_'+str(scaling)+'/', CLARA_dir=CLARA_dir, lattice=self.lattice_file)
        else:
            self.optfunc = partial(self.OptimisingFunction, scaling=scaling, post_injector=self.POST_INJECTOR, CLARA_dir=CLARA_dir, lattice=self.lattice_file)
        # ******************************************************************************
        if not self.POST_INJECTOR:
            best = injector_startingvalues + best
        elif CREATE_BASE_FILES:
            for i in [scaling]:
                pass
                # optfunc(injector_startingvalues + best, scaling=scaling, post_injector=False, verbose=False, runGenesis=False, dir='nelder_mead/basefiles_'+str(i))
        self.set_CLARA_directory(CLARA_dir)
        self.save_parameters = []
        self.changes = None
        # self.base_files = '../../../../basefiles_' + str(int(scaling)) + '/'
        self.verbose = False
        self.best_changes = './transverse_best_changes.yaml'

    def setChangesFile(self, changes):
        self.changes = changes

    def calculateBeamParameters(self):
        # try:
            twiss = self.twiss
            self.framework.change_Lattice_Code('All','elegant')
            self.framework['S02'].prefix = self.base_files
            self.framework.track(startfile='S02', endfile='FMS')

            constraintsList = {}
            quadkls = self.framework.getElementType('quadrupole','k1l')
            quadlengths = self.framework.getElementType('quadrupole','length')
            quadnames = self.framework.getElementType('quadrupole','objectname')
            constraintsListQuads = {
                'max_k': {'type': 'lessthan', 'value': [abs(k/l) for k, l, n in zip(quadkls, quadlengths, quadnames) if not 'FMS' in n], 'limit': 2.5, 'weight': 10},

            }
            constraintsList = merge_two_dicts(constraintsList, constraintsListQuads)

            twiss.read_elegant_twiss_files( [ self.dirname+'/'+n+'.twi' for n in ['S02', 'L02', 'S03', 'L03', 'S04', 'L4H', 'S05', 'S06', 'L04', 'S07', 'FMS']])
            constraintsListSigmas = {
                'max_xrms': {'type': 'lessthan', 'value': 1e3*twiss['sigma_x'], 'limit': 1, 'weight': 5},
                'max_yrms': {'type': 'lessthan', 'value': 1e3*twiss['sigma_y'], 'limit': 1, 'weight': 5},
                'min_xrms': {'type': 'greaterthan', 'value': 1e3*twiss['sigma_x'], 'limit': 0.1, 'weight': 5},
                'min_yrms': {'type': 'greaterthan', 'value': 1e3*twiss['sigma_y'], 'limit': 0.1, 'weight': 5},
                'last_exn': {'type': 'lessthan', 'value': 1e6*twiss['enx'][-1], 'limit': 0.75, 'weight': 100},
                'last_eyn': {'type': 'lessthan', 'value': 1e6*twiss['eny'][-1], 'limit': 0.75, 'weight': 100},
            }
            constraintsList = merge_two_dicts(constraintsList, constraintsListSigmas)

            twiss.read_elegant_twiss_files(self.dirname+'/S07.twi')
            tdc_position = self.framework['CLA-S07-TDC-01-R']['position_start'][2]
            tdc_screen_position = self.framework['CLA-S07-DIA-SCR-03-W']['position_start'][2]
            dechirper_position = self.framework['CLA-S07-DCP-01']['position_start'][2]
            constraintsListS07 = {
                # 'tdc_phase_advance': {'type': 'equalto', 'value': twiss.interpolate(tdc_screen_position,'muy') - twiss.interpolate(tdc_position,'muy'), 'limit': 0.25, 'weight': 0.5},
                # 'tdc_screen_beta_y': {'type': 'greaterthan', 'value': twiss.extract_values('beta_y', tdc_position, tdc_screen_position), 'limit': 5, 'weight': 1},
                'dechirper_sigma_x': {'type': 'lessthan', 'value': 1e3*twiss.interpolate(dechirper_position, 'sigma_x'), 'limit': 0.1, 'weight': 10},
                'dechirper_sigma_y': {'type': 'lessthan', 'value': 1e3*twiss.interpolate(dechirper_position, 'sigma_y'), 'limit': 0.1, 'weight': 10},
                'dechirper_sigma_xy': {'type': 'equalto', 'value': 1e3*twiss.interpolate(dechirper_position, 'sigma_y') - 1e3*twiss.interpolate(dechirper_position, 'sigma_x'), 'limit': 0.0, 'weight': 10},
            }
            constraintsList = merge_two_dicts(constraintsList, constraintsListS07)
            fitness = self.cons.constraints(constraintsList)
            if self.verbose:
                print(self.cons.constraintsList(constraintsList))
            return fitness
        # except:
        #     return 1e6

    def OptimisingFunction(self, inputargs, *args, **kwargs):
        self.inputlist = [a[0]+[a[1]] for a in zip(self.parameters, inputargs)]
        dir = self.optdir+str(self.opt_iteration)
        # print self.inputlist
        fit = self.setup_lattice(self.inputlist, dir, *args, changes=self.changes, **kwargs)
        fitvalue = self.calculateBeamParameters()
        print('fitvalue[', self.opt_iteration, '] = ', fitvalue)
        self.saveState(args, fitvalue)
        self.opt_iteration += 1
        if fitvalue < self.bestfit:
            if hasattr(self, 'constraintsList'):
                print(self.cons.constraintsList(self.constraintsList))
            print('!!!!!!  New best = ', fitvalue)
            pars = (self.parameters+self.save_parameters)
            self.framework.save_changes_file(filename=self.best_changes, elements=pars)
            self.bestfit = fitvalue
        return fitvalue

    def Nelder_Mead(self, best=None, step=0.1):
        best = np.array(self.best) if best is None else np.array(best)
        self.statefile = 'nelder_mead_transverse/best_solutions_running.csv'
        self.optdir = 'nelder_mead_transverse/iteration_'
        print('best = ', best)
        self.bestfit = 1e26

        with open(self.statefile,'w') as out:
            self.opt_iteration = 0
        res = nelder_mead(self.optfunc, best, step=step, max_iter=300, no_improv_break=100)
        print(res)

    def Simplex(self, best=None):
        best = self.best if best is None else best
        self.statefile = 'simplex_transverse/best_solutions_running.csv'
        self.optdir = 'simplex_transverse/iteration_'
        print('best = ', best)
        self.bestfit = 1e26

        with open(self.statefile,'w') as out:
            self.opt_iteration = 0
        res = minimize(self.optfunc, best, method='nelder-mead', options={'disp': True, 'maxiter': 300, 'adaptive': True})
        print(res.x)

    def saveState(self, args, fitness):
        with open(self.statefile,'a') as out:
            csv_out=csv.writer(out)
            args=list(args)
            args.append(self.opt_iteration)
            args.append(fitness)
            csv_out.writerow(args)


if __name__ == "__main__":
        fit = fitnessFunc()
        fit.Nelder_Mead(best)
