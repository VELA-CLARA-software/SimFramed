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

def create_base_files(scaling, charge=70):
    framework = Framework('basefiles_'+str(scaling)+'_'+str(charge), overwrite=True, clean=True)
    framework.loadSettings('CLA10-FE.def')
    if not os.name == 'nt':
        framework.defineASTRACommand(scaling=(scaling))
        framework.defineCSRTrackCommand(scaling=(scaling))
    framework.generator.number_of_particles = 2**(3*(scaling))
    framework.modifyElement('CLA-LRG1-GUN-CAV', 'phase', -5)
    framework.modifyElement('CLA-LRG1-GUN-SOL', 'field_amplitude', 0.15)
    framework['generator'].charge = charge * 1e-12
    # before_tracking()
    framework.track(files=['generator','injector10'])

## Modify as appropriate! ##
# for i in [6]:
    # create_base_files(i,20)
    # create_base_files(i,50)
    # create_base_files(i,70)
    # create_base_files(i,150)
    # create_base_files(i,250)
# exit()

def optFunc(linac1field):
    global dir, lattice
    lattice.modifyElement('CLA-L01-CAV', 'field_amplitude', abs(linac1field[0]))
    lattice.track(startfile='L01', endfile='L01', track=True)
    beam.read_HDF5_beam_file(dir+'/CLA-S02-APER-01.hdf5')
    # print 'mom = ', abs(np.mean(beam.cp) - 38000000)
    return abs(np.mean(beam.cp) - 40500000)


def setMomentum():
    global dir, lattice
    best = [lattice.getElement('CLA-L01-CAV', 'field_amplitude')]
    res = minimize(optFunc, best, method='nelder-mead', tol=1000, options={'disp': False, 'adaptive': True})
    # res = minimize(optFunc, best, method='BFGS', tol=1, options={'eps': 10000, 'disp': True})
    return res.x

def set_Phase(phase, charge=70, track=True):
    global dir, lattice, scaling
    dir = 'SETUP/TOMP_SETUP_'+str(phase)+'_'+str(charge)
    lattice.setSubDirectory(dir)
    # if track:
    #     lattice.change_Lattice_Code('L01','elegant')
    #     lattice.modifyElement('CLA-L01-CAV', 'phase', phase)
    #     lattice['L01'].file_block['input']['prefix'] = '../../../CLARA_BA1_Gaussian/basefiles_'+str(scaling)+'_'+str(charge)+'/'
    #     lattice['L01'].sample_interval = 8
    #     mom = setMomentum()
    #     print('mom = ', mom)
    #     optFunc(mom)
    lattice.change_Lattice_Code('L01','ASTRA')
    lattice['L01'].file_block['input']['prefix'] = '../../../CLARA_BA1_Gaussian/basefiles_'+str(scaling)+'_'+str(charge)+'/'
    lattice['L01'].sample_interval = 8
    lattice.modifyElement('CLA-L01-CAV', 'phase', phase)
    # lattice.modifyElement('CLA-L01-CAV', 'field_amplitude', abs(mom[0]))
    before_tracking()
    if track:
        lattice.track(startfile='L01', endfile='L01')
    lattice['S02'].file_block['input']['prefix'] = '../SETUP/TOMP_SETUP_'+str(phase)+'_'+str(charge)+'/'

def track_phase_elegant(phase, charge=70):
    global dir, lattice
    dir = './NoSC/TOMP_NoSC_'+str(phase)+'_'+str(charge)
    lattice.setSubDirectory(dir)
    lattice.change_Lattice_Code('S02','elegant')
    lattice['S02'].lscDrifts = False
    before_tracking()
    lattice['S02'].file_block['input']['prefix'] = '../../SETUP/TOMP_SETUP_'+str(phase)+'_'+str(charge)+'/'
    lattice.track(startfile='S02', endfile='S02', track=True)

def track_phase_elegant_SC(phase, charge=70):
    global dir, lattice
    dir = './SC/TOMP_SC_'+str(phase)+'_'+str(charge)
    lattice.setSubDirectory(dir)
    lattice.change_Lattice_Code('S02','elegant')
    lattice['S02'].lscDrifts = True
    lattice['S02'].csrDrifts = False
    lattice['S02'].file_block['input']['prefix'] = '../../SETUP/TOMP_SETUP_'+str(phase)+'_'+str(charge)+'/'
    elements = lattice.elementObjects.values()
    for e in elements:
        e.csr_enable = False
        e.sr_enable = False
        e.isr_enable = False
        e.current_bins = 0
        e.lsc_enable = True
        e.lsc_bins = 128
        e.longitudinal_wakefield_enable = True
        e.transverse_wakefield_enable = True
        e.smoothing_half_width = 2
        e.lsc_high_frequency_cutoff_start = 0.25
        e.lsc_high_frequency_cutoff_end = 0.33
    lattices = lattice.latticeObjects.values()
    for l in lattices:
        l.csr_enable = False
        l.lscDrifts = True
        l.csrDrifts = False
        l.lsc_bins = 128
        l.lsc_high_frequency_cutoff_start = 0.25
        l.lsc_high_frequency_cutoff_end = 0.33
    lattice.track(startfile='S02', endfile='S02', track=True)

def track_phase_elegant_SC_CSR(phase, charge=70):
    global dir, lattice
    dir = './SC_CSR/TOMP_SC_CSR_'+str(phase)+'_'+str(charge)
    lattice.setSubDirectory(dir)
    lattice.change_Lattice_Code('S02','elegant')
    lattice['S02'].lscDrifts = True
    lattice['S02'].csrDrifts = False
    lattice['S02'].file_block['input']['prefix'] = '../../SETUP/TOMP_SETUP_'+str(phase)+'_'+str(charge)+'/'
    elements = lattice.elementObjects.values()
    for e in elements:
        e.csr_enable = True
        e.sr_enable = True
        e.isr_enable = True
        e.csr_bins = 50
        e.lsc_enable = True
        e.lsc_bins = 50
        e.longitudinal_wakefield_enable = True
        e.transverse_wakefield_enable = True
        e.smoothing_half_width = 1
        e.lsc_high_frequency_cutoff_start = -1#0.25
        e.lsc_high_frequency_cutoff_end = -1#0.33
    lattices = lattice.latticeObjects.values()
    for l in lattices:
        l.csr_enable = True
        l.lscDrifts = True
        l.csrDrifts = True
        l.lsc_bins = 50
        l.lsc_high_frequency_cutoff_start = -1#0.25
        l.lsc_high_frequency_cutoff_end = -1#0.33
    lattice.track(startfile='S02', endfile='S02', track=True)

def track_phase_astra(phase, charge=70):
    global dir, lattice
    dir = './ASTRA/TOMP_ASTRA_'+str(phase)+'_'+str(charge)
    lattice.setSubDirectory(dir)
    lattice.change_Lattice_Code('S02','ASTRA')
    # before_tracking()
    lattice['S02'].sample_interval = 2**(3*1)
    lattice['S02'].file_block['input']['prefix'] = '../../SETUP/TOMP_SETUP_'+str(phase)+'_'+str(charge)+'/'
    lattice.track(startfile='S02', endfile='S02', track=True)

def track_phase_gpt(phase, charge=70):
    global dir, lattice
    dir = './GPT/TOMP_GPT_'+str(phase)+'_'+str(charge)
    lattice.setSubDirectory(dir)
    lattice.change_Lattice_Code('S02','GPT')
    # before_tracking()
    lattice['S02'].sample_interval = 2**(3*2)
    lattice['S02'].file_block['input']['prefix'] = '../../SETUP/TOMP_SETUP_'+str(phase)+'_'+str(charge)+'/'
    lattices = lattice.latticeObjects.values()
    for l in lattices:
        l.csr_enable = False
        try:
            l.file_block['charge']['space_charge_mode'] = '3D'
        except:
            pass
    lattice.track(startfile='S02', endfile='S02', track=True)

def track_phase_gpt_CSR(phase, charge=70):
    global dir, lattice
    dir = './GPT_CSR/TOMP_GPT_CSR_'+str(phase)+'_'+str(charge)
    lattice.setSubDirectory(dir)
    lattice.change_Lattice_Code('S02','GPT')
    # before_tracking()
    lattice['S02'].sample_interval = 2**(3*1)
    lattice['S02'].file_block['input']['prefix'] = '../../SETUP/TOMP_SETUP_'+str(phase)+'_'+str(charge)+'/'
    lattices = lattice.latticeObjects.values()
    for l in lattices:
        l.csr_enable = True
    lattice.track(startfile='S02', endfile='S02', track=True)

def updateOutput(output):
    sys.stdout.write(output + '\r')
    sys.stdout.flush()

def optFuncVELA(args):
    global dir, lattice, bestdelta
    try:
        lattice.modifyElement('CLA-S02-MAG-QUAD-01', 'k1l', args[0])
        lattice.modifyElement('CLA-S02-MAG-QUAD-02', 'k1l', args[1])
        lattice.modifyElement('CLA-S02-MAG-QUAD-03', 'k1l', args[2])
        lattice.modifyElement('CLA-S02-MAG-QUAD-04', 'k1l', args[3])
        lattice.modifyElement('CLA-S02-MAG-QUAD-05', 'k1l', args[4])
        lattice.track(startfile='S02', endfile='S02', track=True)

        constraintsList = {}
        twiss.reset_dicts()
        twiss.read_elegant_twiss_files(dir+'/S02.twi')
        ipindex = list(twiss['element_name']).index('CLA-VHEE-PHANTOM-01')
        constraintsListVHEE = {
            'ip_betax': {'type': 'lessthan', 'value': twiss['beta_x'][ipindex], 'limit': 0.5, 'weight': 25},
            'ip_betay': {'type': 'lessthan', 'value': twiss['beta_y'][ipindex], 'limit': 0.5, 'weight': 25},
            'ip_betaxy': {'type': 'equalto', 'value': twiss['beta_x'][ipindex], 'limit': twiss['beta_y'][ipindex], 'weight': 25},
            'ip_alphax': {'type': 'equalto', 'value': twiss['alpha_x'][ipindex], 'limit': 0., 'weight': 2.5},
            'ip_alphay': {'type': 'equalto', 'value': twiss['alpha_y'][ipindex], 'limit': 0., 'weight': 2.5},
            'ip_etax': {'type': 'lessthan', 'value': abs(twiss['eta_x'][ipindex]), 'limit': 1e-4, 'weight': 50},
            'ip_etaxp': {'type': 'lessthan', 'value': abs(twiss['eta_xp'][ipindex]), 'limit': 1e-4, 'weight': 50},
        }
        constraintsList = merge_two_dicts(constraintsList, constraintsListVHEE)

        cons = constraintsClass()
        delta = cons.constraints(constraintsList)
        updateOutput('VELA delta = ' + str(delta))
        if delta < bestdelta:
            bestdelta = delta
            print('### New Best: ', delta)
            # print(*args, sep = ", ")
            print '[',', '.join(map(str,args)),']'
            # print(cons.constraintsList(constraintsList))
        # if delta < 0.1:
        #     return 0
    except Exception as e:
        print('Error! = ', e)
        delta = 1e6
    else:
        return delta

def setVELA(quads):
    global dir, lattice, bestdelta
    bestdelta = 1e10
    best = quads
    res = minimize(optFuncVELA, best, method='nelder-mead', options={'disp': False, 'adaptive': True, 'fatol': 1e-3, 'maxiter': 100})
    return res.x

def optimise_Lattice(phase=None, q=70, do_optimisation=False):
    global dir, lattice, scaling, bestdelta
    bestdelta = 1e10
    dir = './SETUP/TOMP_SETUP'
    lattice = Framework(dir, clean=False, verbose=False)
    lattice.loadSettings('CLA10-FE.def')
    lattice.modifyElement('CLA-L01-CAV-SOL-01', 'field_amplitude', 0.07)
    lattice.modifyElement('CLA-L01-CAV-SOL-02', 'field_amplitude', -0.05)
    scaling = 6
    if not os.name == 'nt':
        lattice.defineASTRACommand(scaling=(scaling))
        lattice.defineCSRTrackCommand(scaling=(scaling))
        lattice.define_gpt_command(scaling=(scaling))
    # lattice['L01'].file_block['input']['prefix'] = '../../CLARA_BA1_Gaussian/basefiles_'+str(scaling)+'_'+str(q)+'/'
    quads = [
        lattice.getElement('CLA-S02-MAG-QUAD-01', 'k1l'),
        lattice.getElement('CLA-S02-MAG-QUAD-02', 'k1l'),
        lattice.getElement('CLA-S02-MAG-QUAD-03', 'k1l'),
        lattice.getElement('CLA-S02-MAG-QUAD-04', 'k1l'),
        lattice.getElement('CLA-S02-MAG-QUAD-05', 'k1l'),
    ]
    quads = np.array([  2.018036284043403, -1.6775106326141502, 2.0713066380968943, -1.4092617257511626, 0.07875236608471231 ])
    lattice['S02'].sample_interval = 8
    if do_optimisation:
        if phase is not None:
            lattice['S02'].file_block['input']['prefix'] = '../TOMP_SETUP_'+str(phase)+'_'+str(q)+'/'
        quads = setVELA(quads)
    optFuncVELA(quads)
    print('VELA = ', quads)
    lattice['S02'].sample_interval = 1#2**(3*2)


optimise_Lattice(do_optimisation=False)
# elegantNoSCOut = open('elegant_NoSC_phase_data.csv','w')
# elegantSCOut = open('elegant_SC_phase_data.csv','w')
# elegantSCCSROut = open('elegant_SC_CSR_phase_data.csv','w')
# ASTRAOut = open('ASTRA_phase_data.csv','w')
# GPTout = open('GPT_phase_data.csv','w')
# GPTCSRout = open('GPT_CSR_phase_data.csv','w')

# elegentNoSC_csv_out = csv.writer(elegantNoSCOut)
# elegentSC_csv_out = csv.writer(elegantSCOut)
# elegentSCCSR_csv_out = csv.writer(elegantSCCSROut)
# ASTRA_csv_out = csv.writer(ASTRAOut)
# GPT_csv_out = csv.writer(GPTout)
# GPTCSR_csv_out = csv.writer(GPTCSRout)

for i in [81]:
    for q in [70]:
        # print('q = ', q, '  i = ', i)
        set_Phase(i, q, track=True)
        optimise_Lattice(i, q, True)
        # data = track_phase_elegant(i, q)
        # exit()
        # elegentNoSC_csv_out.writerow(data)
        # data = track_phase_elegant_SC(i, q)
        # elegentSC_csv_out.writerow(data)
        data = track_phase_elegant_SC_CSR(i, q)
        # elegentSCCSR_csv_out.writerow(data)
        # data = track_phase_astra(i, q)
        # # ASTRA_csv_out.writerow(data)
        # data = track_phase_gpt(i)
        # # GPT_csv_out.writerow(data)
        # data = track_phase_gpt_CSR(i)
        # # GPTCSR_csv_out.writerow(data)
