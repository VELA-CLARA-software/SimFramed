import os, errno, sys
import numpy as np
import random
sys.path.append(os.path.abspath(__file__+'/../../../../../'))
from SimulationFramework.Modules.constraints import constraintsClass
from SimulationFramework.Modules.nelder_mead import nelder_mead
import time
import csv
from copy import copy
from functools import partial
from collections import OrderedDict
from shutil import copyfile
from SimulationFramework.Modules.merge_two_dicts import merge_two_dicts
import SimulationFramework.Examples.CLARA.Elegant.runElegant as runEle

def saveState(args, n, fitness):
    with open('nelder_mead/best_solutions_running.csv','a') as out:
        csv_out=csv.writer(out)
        args=list(args)
        args.append(n)
        args.append(fitness)
        csv_out.writerow(args)

def saveParameterFile(best, file='clara_elegant_best.yaml'):
    if POST_INJECTOR:
        allparams = zip(*(parameternames))
    else:
        allparams = zip(*(injparameternames+parameternames))
    output = {}
    for p, k, v in zip(allparams[0], allparams[1], best):
        if p not in output:
            output[p] = {}
        output[p][k] = v
        with open(file,"w") as yaml_file:
            yaml.dump(output, yaml_file, default_flow_style=False)

class Optimise_Elegant(object):

    injector_parameter_names = [
        ['CLA-HRG1-GUN-CAV', 'phase'],
        ['CLA-HRG1-GUN-SOL', 'field_amplitude'],
        ['CLA-L01-CAV', 'field_amplitude'],
        ['CLA-L01-CAV', 'phase'],
        ['CLA-L01-CAV-SOL-01', 'field_amplitude'],
        ['CLA-L01-CAV-SOL-02', 'field_amplitude'],
    ]
    parameter_names = [
        ['CLA-L02-CAV', 'field_amplitude'],
        ['CLA-L02-CAV', 'phase'],
        ['CLA-L03-CAV', 'field_amplitude'],
        ['CLA-L03-CAV', 'phase'],
        ['CLA-L4H-CAV', 'field_amplitude'],
        ['CLA-L4H-CAV', 'phase'],
        ['CLA-L04-CAV', 'field_amplitude'],
        ['CLA-L04-CAV', 'phase'],
        ['bunch_compressor', 'set_angle'],
        ['CLA-S07-DCP-01', 'factor'],
    ]

    def __init__(self):
        super(Optimise_Elegant, self).__init__()
        self.cons = constraintsClass()
        self.changes = None
        self.lattice = None
        self.resultsDict = {}
        # ******************************************************************************
        CLARA_dir = os.path.relpath(__file__+'/../../')
        self.POST_INJECTOR = True
        CREATE_BASE_FILES = False
        self.scaling = 5
        if self.POST_INJECTOR and CREATE_BASE_FILES:
            self.optfunc = partial(self.OptimisingFunction, scaling=self.scaling, post_injector=self.POST_INJECTOR, verbose=False, basefiles='../basefiles_'+str(self.scaling)+'/', CLARA_dir=CLARA_dir)
        else:
            self.optfunc = partial(self.OptimisingFunction, scaling=self.scaling, post_injector=self.POST_INJECTOR, verbose=False, CLARA_dir=CLARA_dir)
        # ******************************************************************************
        if not self.POST_INJECTOR:
            best = injector_startingvalues + best
        elif CREATE_BASE_FILES:
            for i in [self.scaling]:
                pass
                # optfunc(injector_startingvalues + best, scaling=scaling, post_injector=False, verbose=False, runGenesis=False, dir='nelder_mead/basefiles_'+str(i))

    def calculate_constraints(self):
        pass

    def set_changes_file(self, changes):
        self.changes = changes

    def set_lattice_file(self, lattice):
        self.lattice = lattice

    def OptimisingFunction(self, inputargs, *args, **kwargs):
        if not self.POST_INJECTOR:
            parameternames = self.injector_parameter_names + self.parameter_names
        else:
            parameternames = copy(self.parameter_names)
        self.inputlist = map(lambda a: a[0]+[a[1]], zip(parameternames, inputargs))

        self.linac_fields = np.array([i[2] for i in self.inputlist if i[1] == 'field_amplitude'])
        self.linac_phases = np.array([i[2] for i in self.inputlist if i[1] == 'phase'])

        if 'dir' in kwargs.keys():
            dir = kwargs['dir']
            del kwargs['dir']
        else:
            dir = self.optdir+str(self.opt_iteration)

        self.fit = runEle.fitnessFunc(self.lattice)
        self.fit.setup_lattice(self.inputlist, dir, lattice=self.lattice, changes=self.changes, *args, **kwargs)
        fitvalue = self.fit.calculateBeamParameters()
        self.opt_iteration += 1
        constraintsList = self.calculate_constraints()
        fitvalue = self.cons.constraints(constraintsList)
        print 'fitvalue[', self.opt_iteration-1, '] = ', fitvalue
        saveState(inputargs, self.opt_iteration-1, fitvalue)
        if fitvalue < self.bestfit:
            print '!!!!!!  New best = ', fitvalue
            copyfile(dir+'/changes.yaml', self.best_changes)
            self.bestfit = fitvalue
        return fitvalue

    def Nelder_Mead(self, best=None, step=0.1):
        best = np.array(self.best) if best is None else np.array(best)
        self.optdir = 'nelder_mead/iteration_'
        self.best_changes = './nelder_mead_best_changes.yaml'
        print 'best = ', best
        self.bestfit = 1e26

        with open('nelder_mead/best_solutions_running.csv','w') as out:
            self.opt_iteration = 0
        res = nelder_mead(self.optfunc, best, step=step, max_iter=300, no_improv_break=100)
        print res

    def Simplex(self, best=None):
        best = self.best if best is None else best
        self.optdir = 'simplex/iteration_'
        self.best_changes = './simplex_best_changes.yaml'
        print 'best = ', best
        self.bestfit = 1e26

        with open('simplex/best_solutions_running.csv','w') as out:
            self.opt_iteration = 0
        res = minimize(self.optfunc, best, method='nelder-mead', options={'disp': True, 'maxiter': 300, 'adaptive': True})
        print res.x
