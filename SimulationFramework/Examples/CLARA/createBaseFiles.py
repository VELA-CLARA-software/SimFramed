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
        framework.defineASTRACommand(['mpiexec','-np',str(3*scaling),'/opt/ASTRA/astra_MPICH2.sh'])
        framework.defineGeneratorCommand(['/opt/ASTRA/generator.sh'])
        framework.defineCSRTrackCommand(['/opt/OpenMPI-1.4.3/bin/mpiexec','-n',str(3*scaling),'/opt/CSRTrack/csrtrack_openmpi.sh'])
    framework.generator.number_of_particles = 2**(3*scaling)
    framework.track(files=['generator','injector400','S02'])

## Modify as appropriate! ##
for i in [3,4,5,6]:
    create_base_files(i)
exit()
