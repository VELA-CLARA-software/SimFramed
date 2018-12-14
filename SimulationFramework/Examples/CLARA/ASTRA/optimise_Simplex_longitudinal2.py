import numpy as np
import os, sys
sys.path.append(os.path.abspath(__file__+'/../../../../../'))
from FitnessFunc_Longitudinal import *
import operator
import random
import csv
from scipy.optimize import minimize
import yaml
from functools import partial

''' Run the injector part once if only optimising post-injector parameters'''
def create_base_files(scaling):
    framework = Framework('basefiles_'+str(scaling), overwrite=True)
    framework.loadSettings('Lattices/clara400_v12_v3_elegant.def')
    if not os.name == 'nt':
        framework.defineASTRACommand(['mpiexec','-np',str(3*scaling),'/opt/ASTRA/astra_MPICH2.sh'])
        framework.defineGeneratorCommand(['/opt/ASTRA/generator.sh'])
        framework.defineCSRTrackCommand(['/opt/OpenMPI-1.4.3/bin/mpiexec','-n',str(3*scaling),'/opt/CSRTrack/csrtrack_openmpi.sh'])
    framework.generator.number_of_particles = 2**(3*scaling)
    framework.track(files=['S02','L02','S03','L03','S04','L4H','S05','VBC','S06','L04','S07','FMS'])

for i in [3,4]:
    create_base_files(i)
exit()

def saveState(args, fitness):
    with open('best_solutions_running_simplex.csv','a') as out:
        csv_out=csv.writer(out)
        csv_out.writerow(args +[fitness])

def saveParameterFile(best, file='clara_longitudinal_best.yaml'):
    allparams = zip(*(injparameternames+parameternames))
    output = {}
    for p, k, v in zip(allparams[0], allparams[1], best):
        if p not in output:
            output[p] = {}
        output[p][k] = v
        with open(file,"w") as yaml_file:
            yaml.dump(output, yaml_file, default_flow_style=False)

def optfunc(inputargs, dir=None, *args, **kwargs):
    if dir == None:
        with TemporaryDirectory(dir=os.getcwd()) as tmpdir:
            fit = fitnessFunc(inputargs, tmpdir, *args, **kwargs)
            fitvalue = fit.calculateBeamParameters()
    else:
            fit = fitnessFunc(inputargs, dir, *args, **kwargs)
            fitvalue = fit.calculateBeamParameters()
    saveState(inputargs, fitvalue)
    return fitvalue

framework = Framework('longitudinal_best', overwrite=False)
framework.loadSettings('Lattices/clara400_v12_v3.def')
injparameters = []
parameters = []
injparameternames = [['CLA-HRG1-GUN-CAV', 'phase'], ['CLA-HRG1-GUN-SOL', 'field_amplitude'], ['CLA-L01-CAV', 'field_amplitude'],
                     ['CLA-L01-CAV', 'phase'], ['CLA-L01-CAV-SOL-01', 'field_amplitude'], ['CLA-L01-CAV-SOL-02', 'field_amplitude']]
parameternames = [
['CLA-L02-CAV', 'field_amplitude'], ['CLA-L02-CAV', 'phase'], ['CLA-L03-CAV', 'field_amplitude'], ['CLA-L03-CAV', 'phase'],
['CLA-L4H-CAV', 'field_amplitude'], ['CLA-L4H-CAV', 'phase'], ['CLA-L04-CAV', 'field_amplitude'], ['CLA-L04-CAV', 'phase'],
['bunch_compressor','angle']
]
for p in injparameternames:
    injparameters.append(framework.getElement(*p))

''' always '''
for p in parameternames:
    parameters.append(framework.getElement(*p))

results = []
with open('CLARA_longitudinal_best_solutions_simplex.csv.tmp', 'r') as csvfile:
  reader = csv.reader(csvfile, delimiter=',', quoting=csv.QUOTE_NONNUMERIC)
  for row in reader:
    results.append(row)
best = results[0]
print 'starting values = ', best

''' export best results to a yaml file '''
# saveParameterFile(best)
# exit()
# print 'best fitness = ', optfunc(best, dir=os.getcwd()+'/CLARA_best_longitudinal_simplex', scaling=6, overwrite=True, verbose=True, summary=True, post_injector=False)
# exit()

''' if including injector'''
# best = injparameters + parameters
''' ELSE '''
# best = parameters

optfunc3 = partial(optfunc, dir=None, scaling=3, post_injector=False)
optfunc4 = partial(optfunc, dir=None, scaling=4, post_injector=False)

with open('best_solutions_running_simplex.csv','w') as out:
    pass
res = minimize(optfunc3, best, method='nelder-mead', options={'xtol': 1e-8, 'disp': True, 'adaptive': True})
print res

try:
    print 'best fitness = ', optfunc(res.x, dir=os.getcwd()+'/CLARA_best_longitudinal_simplex', scaling=6, overwrite=True, verbose=True, summary=True, post_injector=False)
except:
    pass

saveParameterFile(res.x)
