import numpy as np
import os, sys
sys.path.append(os.path.abspath(__file__+'/../../../../'))
from SimulationFramework.Framework import *
import csv
import shutil
import uuid

''' Run the injector part once if only optimising post-injector parameters'''
def create_base_files(scaling):
    framework = Framework('basefiles_'+str(scaling), overwrite=True)
    framework.loadSettings('Lattices/clara400_v12_v3.def')
    if not os.name == 'nt':
        framework.defineASTRACommand(scaling=scaling)
        framework.defineCSRTrackCommand(scaling=scaling)
    framework.generator.number_of_particles = 2**(3*scaling)
    # framework.change_Lattice_Code('VBC', 'ASTRA')
    framework.track()#startfile='VBC')

## Modify as appropriate! ##
for i in [4,5,6]:
    create_base_files(i)
exit()
