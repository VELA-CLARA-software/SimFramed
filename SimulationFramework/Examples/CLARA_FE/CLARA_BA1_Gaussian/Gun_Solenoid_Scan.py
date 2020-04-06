import sys, os
import numpy as np
import csv
from scipy.optimize import minimize
sys.path.append('../../../../')
from SimulationFramework.Framework import *
import SimulationFramework.Modules.read_beam_file as rbf
beam = rbf.beam()
import SimulationFramework.Modules.read_twiss_file as rtf
from SimulationFramework.Modules.merge_two_dicts import merge_two_dicts
from SimulationFramework.Modules.constraints import *
twiss = rtf.twiss()


def before_tracking():
    elements = lattice.elementObjects.values()
    for e in elements:
        e.lsc_enable = True
        e.lsc_bins = 200
        e.current_bins = 0
        e.longitudinal_wakefield_enable = True
        e.transverse_wakefield_enable = True
    lattices = lattice.latticeObjects.values()
    for l in lattices:
        l.lscDrifts = True
        l.lsc_bins = 200

def sol_field_from_current(I):
    return 0.001203*I + 1.09776e-7*I**2 - 2.53708e-10*I**3

def bsol_field_from_current(I):
    return 0.354*sol_field_from_current(I)

def create_base_files(scaling, charge=70, sol=160, bsol=-144, gunfield=58.82):
    sol=int(sol)
    bsol=int(bsol)
    framework = Framework('Gun_Solenoid_Scan/basefiles_'+str(scaling)+'_'+str(charge)+'_'+str(int(sol))+'_'+str(int(bsol))+'_'+str(int(gunfield))+'MVm', overwrite=True, clean=True, verbose=False)
    framework.loadSettings('CLA10-BA1_TOMP.def')
    if not os.name == 'nt':
        framework.defineASTRACommand(scaling=(scaling))
        framework.defineCSRTrackCommand(scaling=(scaling))
    framework.generator.number_of_particles = 2**(3*scaling)
    framework.modifyElement('CLA-LRG1-GUN-CAV', 'phase', 0)
    framework.modifyElement('CLA-LRG1-GUN-CAV', 'field_amplitude', gunfield*1e6) #66000000 == 4.7MeV/c   58820000 == 4.3MeV/c
    framework.modifyElement('CLA-LRG1-MAG-SOL-01', 'field_amplitude', sol_field_from_current(sol))
    framework.modifyElement('CLA-LRG1-MAG-BSOL-01', 'field_amplitude', bsol_field_from_current(bsol))
    framework.modifyElement('CLA-L01-CAV', 'phase', 0)
    framework.modifyElement('CLA-L01-CAV', 'field_amplitude', 21200000)
    framework.modifyElement('CLA-L01-MAG-SOL-01', 'field_amplitude', 0)
    framework.modifyElement('CLA-L01-MAG-SOL-02', 'field_amplitude', 0)
    framework.modifyElement('CLA-S02-MAG-QUAD-01', 'k1l', 0)
    framework.modifyElement('CLA-S02-MAG-QUAD-02', 'k1l', 0)
    framework.modifyElement('CLA-S02-MAG-QUAD-03', 'k1l', 0)
    framework.modifyElement('CLA-S02-MAG-QUAD-04', 'k1l', 0)
    framework.generator.load_defaults('clara_400_2ps_Gaussian')
    framework['generator'].charge = charge * 1e-12
    framework.change_Lattice_Code('L01','ASTRA')
    framework.change_Lattice_Code('S02','ASTRA')
    framework.track(files=['generator', 'injector10', 'L01', 'S02'])

# Modify as appropriate! ##
for i in [4]:
    print('Scaling = ', i)
    for j in np.arange(120,251,1):
        print('  Solenoid = ', j)
        for k in [53,55,57]:#np.arange(50,70,2):
            print('    Gun Field = ', k)
            create_base_files(i, 70, j, -0.9*j, gunfield=k)
    # for j in np.arange(181,220,1):
    #     print('  Solenoid = ', j)
    #     create_base_files(i, 70, j, -0.9*j)
# for i in [4]:
#     print('Scaling = ', i)
#     for j in np.arange(120,121,10):
#         print('  Solenoid = ', j)
#         create_base_files(i, 70, j, -0.9*j)
exit()
