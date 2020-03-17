import numpy as np
import os, sys
sys.path.append(os.path.abspath(__file__+'/../../../../'))
from SimulationFramework.Framework import *
import csv
import shutil
import uuid
import SimulationFramework.Modules.read_twiss_file as rtf
twiss = rtf.twiss()

''' Run the injector part once if only optimising post-injector parameters'''
def create_base_files(scaling):
    framework = Framework('basefiles_'+str(scaling), overwrite=True)
    framework.loadSettings('Lattices/clara400_v12_v3.def')
    if not os.name == 'nt':
        framework.defineASTRACommand(scaling=scaling)
        framework.defineCSRTrackCommand(scaling=scaling)
    framework.generator.number_of_particles = 2**(3*scaling)
    # framework.change_Lattice_Code('VBC', 'ASTRA')
    framework.track(endfile="S07")#startfile='VBC')

''' Run the injector part once if only optimising post-injector parameters'''
def create_base_files_20pC(scaling):
    framework = Framework('basefiles_'+str(scaling)+'_20pC', overwrite=True)
    framework.loadSettings('Lattices/clara400_v12_v3.def')
    if not os.name == 'nt':
        framework.defineASTRACommand(scaling=scaling)
        framework.defineCSRTrackCommand(scaling=scaling)
    framework.generator.number_of_particles = 2**(3*scaling)
    framework.generator.charge = 20e-12
    framework['CLA-HRG1-GUN-SOL'].field_amplitude = 0.322
    framework['CLA-L01-CAV-SOL-01'].field_amplitude = 0.08
    framework['CLA-L01-CAV-SOL-02'].field_amplitude = -0.08
    framework['CLA-L02-CAV'].phase = 20
    framework['CLA-L03-CAV'].phase = 30
    # framework['bunch_compressor'].angle = 0.1
    # framework['L02'].sample_interval = 2**(3*2)
    framework.track(endfile="S07")#, startfile='VBC')

''' Run the injector part once if only optimising post-injector parameters'''
def create_base_files_100pC(scaling, sol=0.34, lsol=0.08, subdir=True):
    if subdir:
        dir = 'basefiles_'+str(scaling)+'_100pC/'+str(sol)+'/'+str(lsol)+'/'
    else:
        dir = 'basefiles_'+str(scaling)+'_100pC/'
    framework = Framework(dir, overwrite=True, clean=True)
    if subdir:
        framework.loadSettings('Lattices/clara400_v12_v3_elegant.def')
    else:
        framework.loadSettings('Lattices/clara400_v12_v3.def')
    if not os.name == 'nt':
        framework.defineASTRACommand(scaling=scaling)
        framework.defineCSRTrackCommand(scaling=scaling)
    framework.generator.number_of_particles = 2**(3*scaling)
    framework.generator.charge = 100e-12
    framework['CLA-HRG1-GUN-SOL'].field_amplitude = sol
    framework['CLA-L01-CAV-SOL-01'].field_amplitude = lsol
    framework['CLA-L01-CAV-SOL-02'].field_amplitude = -lsol
    framework['CLA-L02-CAV'].phase = -15
    framework['CLA-L03-CAV'].phase = -12
    framework['bunch_compressor'].angle = 0.135
    # framework['L02'].sample_interval = 2**(3*2)
    framework.track(endfile="S07")#, startfile='VBC')
    twiss.read_astra_emit_files(dir + '/injector400.Xemit.001')
    return twiss['enx'][-1]

## Modify as appropriate! ##
# for i in [4]:
#     for j in np.arange(0.321,0.331,0.001):
#         for k in np.arange(0,0.11,0.02):
#             enx = create_base_files_100pC(i,j,k)
#             print('SOL = ', j, '  BSOL = ', k, '  enx = ', enx)
for i in [4,5,6]:
    enx = create_base_files_100pC(scaling=i,sol=0.323,lsol=0.0, subdir=False)
exit()
