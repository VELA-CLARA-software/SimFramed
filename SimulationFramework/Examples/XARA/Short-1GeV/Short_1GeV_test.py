import sys, os
import numpy as np
from Short_1GeV_nelder_mead import Short1GeV

if __name__ == "__main__":

    opt = Short1GeV()
    opt.set_changes_file(['../GA_best_changes.yaml', '../transverse_best_changes.yaml'])
    opt.set_lattice_file('Lattices/claraX400_v12_80MVm_Elegant.def')
    opt.base_files = '../../../CLARA/basefiles_6/'
    opt.start_lattice = 'CLARAX'
    opt.optdir = 'test'
    opt.opt_iteration = ''
    opt.runElegant = True
    opt.runGenesis = False
    opt.best_changes = './test_best_changes.yaml'
    opt.sample_interval=2**(3*1)
    opt.ncpu = 24
    opt.elegant_ncpu = 1
    opt.genesis_ncpu = 22
    opt.delete_most_output = False
    best = opt.load_best('./GA_best_changes.yaml')
    opt.OptimisingFunction(best)
