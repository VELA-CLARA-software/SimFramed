from ASTRAInjector import *
import numpy as np
from constraints import *
import os
import tempfile
import copy
import read_twiss_file as rtf
from deap import algorithms
from deap import base
from deap import creator
from deap import tools
import multiprocessing
from scoop import futures
import operator
import random

import csv

def merge_two_dicts(x, y):
    """Given two dicts, merge them into a new dict as a shallow copy."""
    z = x.copy()
    z.update(y)
    return z


twiss = rtf.twiss()

twiss.read_astra_emit_files(['twiss_best_2120/test.1.Zemit.001', 'twiss_best_2120/test.2.Zemit.001',
                             'twiss_best_2120/test.3.Zemit.001', 'twiss_best_2120/test.4.Zemit.001',
                             'twiss_best_2120/test.5.Zemit.001'])
# print twiss['mux']
print twiss.interpolate(46.4378,'muy') - twiss.interpolate(40.8012,'muy')
# print twiss.extract_values('beta_x', 40.8012, 46.4378)
print twiss.interpolate(46.4378, 'beta_x')
dirname = 'twiss_best_2120'
cons = constraintsClass()
constraintsList = {}
twiss.read_astra_emit_files(dirname+'/test.2.Zemit.001')
constraintsList2 = {
    'max_xrms_2': {'type': 'lessthan', 'value': 1e3*max(twiss['sigma_x']), 'limit': 1, 'weight': 10},
    'max_yrms_2': {'type': 'lessthan', 'value': 1e3*max(twiss['sigma_y']), 'limit': 1, 'weight': 10},
    'min_xrms_2': {'type': 'greaterthan', 'value': 1e3*min(twiss['sigma_x']), 'limit': 0.2, 'weight': 10},
    'min_yrms_2': {'type': 'greaterthan', 'value': 1e3*min(twiss['sigma_y']), 'limit': 0.2, 'weight': 10},
    'beta_x_2': {'type': 'lessthan', 'value': twiss['beta_x'], 'limit': 50, 'weight': 50},
    'beta_y_2': {'type': 'lessthan', 'value': twiss['beta_y'], 'limit': 50, 'weight': 50},
}
constraintsList = merge_two_dicts(constraintsList, constraintsList2)
# print constraintsList2
twiss.read_astra_emit_files(dirname+'/test.3.Zemit.001')
constraintsList3 = {
    'max_xrms_3': {'type': 'lessthan', 'value': 1e3*max(twiss['sigma_x']), 'limit': 1, 'weight': 10},
    'max_yrms_3': {'type': 'lessthan', 'value': 1e3*max(twiss['sigma_y']), 'limit': 1, 'weight': 10},
    'min_xrms_3': {'type': 'greaterthan', 'value': 1e3*min(twiss['sigma_x']), 'limit': 0.2, 'weight': 10},
    'min_yrms_3': {'type': 'greaterthan', 'value': 1e3*min(twiss['sigma_y']), 'limit': 0.2, 'weight': 10},
    'beta_x_3': {'type': 'lessthan', 'value': twiss['beta_x'], 'limit': 50, 'weight': 50},
    'beta_y_3': {'type': 'lessthan', 'value': twiss['beta_y'], 'limit': 50, 'weight': 50},
}
constraintsList = merge_two_dicts(constraintsList, constraintsList3)
# print constraintsList3
twiss.read_astra_emit_files(dirname+'/test.4.Zemit.001')
constraintsList4 = {
    'min_xrms_4': {'type': 'greaterthan', 'value': 1e3*min(twiss['sigma_x']), 'limit': 0.1, 'weight': 30},
    'min_yrms_4': {'type': 'greaterthan', 'value': 1e3*min(twiss['sigma_y']), 'limit': 0.2, 'weight': 50},
    'last_yrms_4': {'type': 'greaterthan', 'value': 1e3*twiss['sigma_y'][-1], 'limit': 0.4, 'weight': 30},
    'last_xrms_4': {'type': 'lessthan', 'value': 1e3*twiss['sigma_x'][-1], 'limit': 0.2, 'weight': 30},
    'last_exn_4': {'type': 'lessthan', 'value': 1e6*twiss['enx'][-1], 'limit': 0.8, 'weight': 10},
    'last_eyn_4': {'type': 'lessthan', 'value': 1e6*twiss['eny'][-1], 'limit': 0.8, 'weight': 10},
    'beta_x_4': {'type': 'lessthan', 'value': twiss['beta_x'], 'limit': 50, 'weight': 50},
    'beta_y_4': {'type': 'lessthan', 'value': twiss['beta_y'], 'limit': 50, 'weight': 50},
}
constraintsList = merge_two_dicts(constraintsList, constraintsList4)
# print constraintsList4
twiss.read_astra_emit_files(dirname+'/test.5.Zemit.001')
constraintsList5 = {
    'min_xrms_5': {'type': 'greaterthan', 'value': 1e3*min(twiss['sigma_x']), 'limit': 0.1, 'weight': 20},
    'min_yrms_5': {'type': 'greaterthan', 'value': 1e3*min(twiss['sigma_y']), 'limit': 0.2, 'weight': 20},
    'last_alpha_x_5': {'type': 'lessthan', 'value': abs(twiss['alpha_x'][-1]), 'limit': 2, 'weight': 10},
    'last_alpha_y_5': {'type': 'lessthan', 'value': abs(twiss['alpha_y'][-1]), 'limit': 2, 'weight': 10},
    'last_beta_x_5': {'type': 'lessthan', 'value': twiss['beta_x'], 'limit': 50, 'weight': 50},
    'tdc_phase_advance': {'type': 'equalto', 'value': twiss.interpolate(46.4378,'muy') - twiss.interpolate(40.8012,'muy'), 'limit': 0.25, 'weight': 500},
    'tdc_beta_y_greaterthan': {'type': 'greaterthan', 'value': twiss.interpolate(40.8012, 'beta_y'), 'limit': 50, 'weight': 25},
    'tdc_beta_y_lassthan': {'type': 'lessthan', 'value': twiss.interpolate(40.8012, 'beta_y'), 'limit': 100, 'weight': 25},
    'tdc_screen_beta_y': {'type': 'greaterthan', 'value': twiss.extract_values('beta_y', 40.8012, 46.4378), 'limit': 3, 'weight': 50},
    'screen_beta_x': {'type': 'equalto', 'value': twiss.interpolate(46.4378, 'beta_x'), 'limit': 5, 'weight': 50},
    'screen_beta_y': {'type': 'equalto', 'value': twiss.interpolate(46.4378, 'beta_y'), 'limit': 5, 'weight': 50},
}
constraintsList = merge_two_dicts(constraintsList, constraintsList5)
# print constraintsList5
print twiss.interpolate(46.4378, 'beta_y')
print cons.constraints(constraintsList5)
