import numpy as np
import os, sys
sys.path.append('../../../../')
import SimulationFramework.Framework as fw
import csv
import shutil
import uuid

''' Run the injector part once if only optimising post-injector parameters'''
def create_base_files(scaling):
    framework = fw.Framework('basefiles_'+str(scaling), clean=False)
    framework.loadSettings('CLA10-BA1_TOMP.def')
    # framework.loadSettings('Lattices/clara10_v12.def')
    if not os.name == 'nt':
        framework.defineASTRACommand(['mpiexec','-np',str(4*scaling),'/opt/ASTRA/astra_MPICH2.sh'])
        framework.defineGeneratorCommand(['/opt/ASTRA/generator.sh'])
        framework.defineCSRTrackCommand(['/opt/OpenMPI-1.4.3/bin/mpiexec','-n',str(4*scaling),'/opt/CSRTrack/csrtrack_openmpi.sh'])
    framework.defineElegantCommand(['elegant'])
    framework.generator.number_of_particles = 2**(3*scaling)
    framework.modifyElement('CLA-LRG1-GUN-CAV', 'phase', -5)
    # framework.modifyElement('CLA-HRG1-GUN-SOL', 'field_amplitude', gunsol)
    # framework.modifyElement('CLA-L01-CAV', 'field_amplitude', abs(linac1field))
    framework.modifyElement('CLA-L01-CAV', 'phase', +10)
    # framework.track(files=['generator','injector10'])#startfile='VBC')
    framework.track(files=['L01','S02','C2V','VELA'])

## Modify as appropriate! ##
for i in [4]:
    create_base_files(i)
exit()
