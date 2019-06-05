import sys, os
import numpy as np
sys.path.append('./../')
from Optimise_Genesis_Elegant import Optimise_Genesis_Elegant
from Short_1GeV_nelder_mead import Short1GeV
from functools import partial

if __name__ == "__main__":

    opt = Short1GeV()
    opt.set_changes_file(['../nelder_mead_best_changes.yaml', '../transverse_best_changes.yaml'])
    opt.set_lattice_file('Lattices/claraX400_v12_80MVm_Elegant.def')
    opt.base_files = '../../../CLARA/basefiles_6/'
    opt.start_lattice = 'CLARAX'
    opt.optdir = 'test'
    opt.opt_iteration = ''
    opt.runElegant = True
    opt.runGenesis = True
    opt.best_changes = './test_best_changes.yaml'
    best = opt.load_best('../nelder_mead_best_changes.yaml')
    opt.OptimisingFunction(best)
