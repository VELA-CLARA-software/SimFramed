import numpy as np
import os, sys
sys.path.append(os.path.abspath(__file__+'/../../../../'))
from SimulationFramework.Framework import *
import FitnessFunc_Longitudinal as FF
from deap import algorithms
from deap import base
from deap import creator
from deap import tools
import operator
import random
import multiprocessing
from scoop import futures
import csv
from SimulationFramework.Modules.optimiser import optimiser
opt = optimiser()
import shutil
import uuid

''' Run the injector part once if only optimising post-injector parameters'''
def create_base_files(scaling):
    framework = Framework('basefiles_'+str(scaling), overwrite=True)
    framework.loadSettings('Lattices/clara400_v12_v3.def')
    if not os.name == 'nt':
        framework.defineASTRACommand(['mpiexec','-np',str(3*scaling),'/opt/ASTRA/astra_MPICH2.sh'])
        framework.defineGeneratorCommand(['/opt/ASTRA/generator.sh'])
        framework.defineCSRTrackCommand(['/opt/OpenMPI-1.4.3/bin/mpiexec','-n',str(3*scaling),'/opt/CSRTrack/csrtrack_openmpi.sh'])
    framework.generator.number_of_particles = 2**(3*scaling)
    framework.track(files=['generator','injector400','S02'])

# for i in [3,4,5,6]:
#     create_base_files(i)
# exit()


def optfunc(args, dir=None, **kwargs):
    if dir == None:
        with FF.TemporaryDirectory(dir=os.getcwd()) as tmpdir:
            fit = FF.fitnessFunc(args, tmpdir, **kwargs)
            fitvalue = fit.calculateBeamParameters()
    else:
            fit = FF.fitnessFunc(args, dir, **kwargs)
            fitvalue = fit.calculateBeamParameters()
    return (fitvalue,)

framework = Framework('longitudinal_best', overwrite=False)
framework.loadSettings('Lattices/clara400_v12_v3.def')
parameters = []
''' if including injector'''
# parameters.append(framework.getElement('CLA-HRG1-GUN-CAV', 'phase'))
# parameters.append(framework.getElement('CLA-HRG1-GUN-SOL', 'field_amplitude'))
# parameters.append(framework.getElement('CLA-L01-CAV', 'field_amplitude'))
# parameters.append(framework.getElement('CLA-L01-CAV', 'phase'))
# parameters.append(framework.getElement('CLA-L01-CAV-SOL-01', 'field_amplitude'))
# parameters.append(framework.getElement('CLA-L01-CAV-SOL-02', 'field_amplitude'))
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
# parameters.append(0.15)
best = parameters

results = []
with open('CLARA_longitudinal_best_solutions_simplex.csv.tmp', 'r') as csvfile:
  reader = csv.reader(csvfile, delimiter=',', quoting=csv.QUOTE_NONNUMERIC)
  for row in reader:
    results.append(row)
best = results[0][-len(parameters):]

best = [x  -7.07101854925 for x in [22200198.49138267,-14.219929641344065,21501627.591965452,0.45153757793072558,23141833.859653812,161.16443024124007,29079460.798801307,39.477777098720388,7.2856647892383979]]

print 'starting values = ', best
# exit()
fit = FF.fitnessFunc(best, os.getcwd()+'/test_6', scaling=6, overwrite=True, verbose=True, summary=False)
print fit.calculateBeamParameters()
exit()
