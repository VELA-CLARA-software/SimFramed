import sys, os
import numpy as np
import csv
from scipy.optimize import minimize
sys.path.append('../../../')
from SimulationFramework.Framework import *
import SimulationFramework.Modules.read_beam_file as rbf
beam = rbf.beam()
import SimulationFramework.Modules.read_twiss_file as rtf
twiss = rtf.twiss()
sys.path.append('../../../../')
import Software.Utils.unitsConverter.unitsConverter as unitsC

def create_base_files(scaling):
    framework = Framework('basefiles_'+str(scaling), overwrite=True, clean=True)
    framework.loadSettings('CLA10-BA1_TOMP.def')
    if not os.name == 'nt':
        framework.defineASTRACommand(['mpiexec','-np',str(3*scaling),'/opt/ASTRA/astra_MPICH2.sh'])
        framework.defineGeneratorCommand(['/opt/ASTRA/generator.sh'])
        framework.defineCSRTrackCommand(['/opt/OpenMPI-1.4.3/bin/mpiexec','-n',str(3*scaling),'/opt/CSRTrack/csrtrack_openmpi.sh'])
    framework.generator.number_of_particles = 2**(3*scaling)
    framework.modifyElement('CLA-LRG1-GUN-CAV', 'phase', -5)
    framework.track(files=['generator','injector10'])

# Modify as appropriate! ##
# for i in [4]:
#     create_base_files(i)
# exit()

def optFunc(linac1field):
    global dir, lattice
    lattice.modifyElement('CLA-L01-CAV', 'field_amplitude', abs(linac1field[0]))
    lattice.track(startfile='L01', endfile='L01', track=True)
    time.sleep(0.1)
    beam.read_HDF5_beam_file(dir+'/CLA-S02-APER-01.hdf5')
    # print 'mom = ', abs(np.mean(beam.cp) - 38000000)
    return abs(np.mean(beam.cp) - 35500000)


def setMomentum():
    global dir, lattice
    best = [lattice.getElement('CLA-L01-CAV', 'field_amplitude')]
    res = minimize(optFunc, best, method='nelder-mead', tol=1000, options={'disp': False, 'adaptive': True})
    # res = minimize(optFunc, best, method='BFGS', tol=1, options={'eps': 10000, 'disp': True})
    return res.x

def set_Phase(phase):
    global dir, lattice
    dir = 'TOMP_SETUP_'+str(phase)
    # lattice.setSubDirectory(dir)
    # lattice.change_Lattice_Code('L01','elegant')
    # lattice.modifyElement('CLA-L01-CAV', 'phase', -1*phase)
    # mom = setMomentum()
    # optFunc(mom)
    # lattice.change_Lattice_Code('L01','ASTRA')
    # lattice.track(startfile='L01', endfile='L01', track=True)
    lattice['S02'].file_block['input']['prefix'] = '../TOMP_SETUP_'+str(phase)+'/'

def track_phase_elegant(phase):
    global dir, lattice
    dir = './TOMP_NoSC_DBURT_'+str(phase)
    lattice.setSubDirectory(dir)
    lattice.change_Lattice_Code('S02','elegant')
    lattice['S02'].lscDrifts = False
    lattice.change_Lattice_Code('C2V','elegant')
    lattice['C2V'].lscDrifts = False
    lattice.change_Lattice_Code('VELA','elegant')
    lattice['VELA'].lscDrifts = False
    lattice.track(startfile='S02', endfile='S02', track=True)
    beam.read_HDF5_beam_file(dir+'/EBT-INJ-DIA-YAG-05.hdf5')
    ansNear = [phase, np.mean(beam.cp), beam.sigma_z, beam.momentum_spread]
    print 'Phase = ', phase
    print '\t\tElegant No SC NEAR:  sigma_z = ', beam.sigma_z, '    momentum-spread = ', beam.momentum_spread
    beam.read_HDF5_beam_file(dir+'/EBT-BA1-COFFIN-FOC.hdf5')
    ansFar = [beam.sigma_z, beam.momentum_spread]
    print '\t\tElegant No SC FAR:   sigma_z = ', beam.sigma_z, '    momentum-spread = ', beam.momentum_spread
    ans = ansNear + ansFar
    return ans

def track_phase_elegant_SC(phase):
    global dir, lattice
    dir = './TOMP_SC_DBURT_'+str(phase)
    lattice.setSubDirectory(dir)
    lattice.change_Lattice_Code('S02','elegant')
    lattice['S02'].lscDrifts = True
    lattice.change_Lattice_Code('C2V','elegant')
    lattice['C2V'].lscDrifts = True
    lattice.change_Lattice_Code('VELA','elegant')
    lattice['VELA'].lscDrifts = True
    lattice.track(startfile='S02', endfile='S02', track=True)
    beam.read_HDF5_beam_file(dir+'/EBT-INJ-DIA-YAG-05.hdf5')
    ansNear = [phase, np.mean(beam.cp), beam.sigma_z, beam.momentum_spread]
    print 'Phase = ', phase
    print '\t\tElegant SC NEAR:  sigma_z = ', beam.sigma_z, '    momentum-spread = ', beam.momentum_spread
    beam.read_HDF5_beam_file(dir+'/EBT-BA1-COFFIN-FOC.hdf5')
    ansFar = [beam.sigma_z, beam.momentum_spread]
    print '\t\tElegant SC FAR:   sigma_z = ', beam.sigma_z, '    momentum-spread = ', beam.momentum_spread
    ans = ansNear + ansFar
    return ans

def track_phase_astra(phase):
    global dir, lattice
    dir = './TOMP_ASTRA_DBURT_'+str(phase)
    lattice.setSubDirectory(dir)
    lattice.change_Lattice_Code('S02','ASTRA')
    lattice.track(startfile='S02', endfile='S02', track=True)
    beam.read_HDF5_beam_file(dir+'/EBT-INJ-DIA-YAG-05.hdf5')
    ansNear = [phase, np.mean(beam.cp), beam.sigma_z, beam.momentum_spread]
    print 'Phase = ', phase
    print '\t\tASTRA NEAR:  sigma_z = ', beam.sigma_z, '    momentum-spread = ', beam.momentum_spread
    beam.read_HDF5_beam_file(dir+'/EBT-BA1-COFFIN-FOC.hdf5')
    ansFar = [beam.sigma_z, beam.momentum_spread]
    print '\t\tASTRA FAR:   sigma_z = ', beam.sigma_z, '    momentum-spread = ', beam.momentum_spread
    ans = ansNear + ansFar
    return ans

def updateOutput(output):
    sys.stdout.write(output + '\r')
    sys.stdout.flush()

def optFuncChicane(args):
    global dir, lattice
    if len(args) == 6:
        s02q1, s02q2, s02q3, s02q4, c2vq1, c2vq2 = args
        lattice.modifyElement('CLA-C2V-MAG-QUAD-01', 'k1', c2vq1)
        lattice.modifyElement('CLA-C2V-MAG-QUAD-02', 'k1', c2vq2)
        lattice.modifyElement('CLA-C2V-MAG-QUAD-03', 'k1', c2vq1)
    else:
        s02q1, s02q2, s02q3, s02q4= args
    lattice.modifyElement('CLA-S02-MAG-QUAD-01', 'k1', s02q1)
    lattice.modifyElement('CLA-S02-MAG-QUAD-02', 'k1', s02q2)
    lattice.modifyElement('CLA-S02-MAG-QUAD-03', 'k1', s02q3)
    lattice.modifyElement('CLA-S02-MAG-QUAD-04', 'k1', s02q4)

    lattice.track(startfile='S02', endfile='C2V', track=True)
    beam.read_HDF5_beam_file(dir+'/CLA-C2V-MARK-01.hdf5')
    twiss.read_sdds_file(dir+'/C2V.twi')
    emitStart = beam.normalized_horizontal_emittance
    betaXStart = twiss['betax'][0]
    betaYStart = twiss['betay'][0]
    # print 'Start beta_x = ', betaXStart, '  beta_y = ', betaYStart, '  emit_x = ', emitStart
    beam.read_HDF5_beam_file(dir+'/CLA-C2V-MARK-02.hdf5')
    emitEnd = beam.normalized_horizontal_emittance
    betaXEnd = twiss['betax'][-1]
    betaYEnd = twiss['betay'][-1]
    betaYPenalty = max(twiss['betay']) - 50 if max(twiss['betay']) > 50 else 0
    etaXEnd = twiss['etax'][-1]
    etaXEnd = etaXEnd if abs(etaXEnd) > 1e-3 else 0
    etaXPEnd = twiss['etaxp'][-1]
    etaXPEnd = etaXEnd if abs(etaXEnd) > 1e-3 else 0
    # print 'End beta_x = ', betaXEnd, '  beta_y = ', betaYEnd, '  emit_x = ', emitEnd
    delta = np.sqrt((1e6*abs(emitStart-emitEnd))**2 + 1*abs(betaXStart-betaXEnd)**2 + 1*abs(betaYStart-betaYEnd)**2 + 10*abs(betaYPenalty)**2 + 100*abs(etaXEnd)**2 + 100*abs(etaXPEnd)**2)
    if delta < 0.1:
        delta = 0.1
    updateOutput('Chicane delta = ' + str(delta))
    return delta

def setChicane(quads):
    global dir, lattice
    # best = [lattice.getElement('CLA-S02-MAG-QUAD-01', 'k1'),
    #         lattice.getElement('CLA-S02-MAG-QUAD-02', 'k1'),
    #         lattice.getElement('CLA-S02-MAG-QUAD-03', 'k1'),
    #         lattice.getElement('CLA-S02-MAG-QUAD-04', 'k1'),
    #         lattice.getElement('CLA-C2V-MAG-QUAD-01', 'k1'),
    #         lattice.getElement('CLA-C2V-MAG-QUAD-02', 'k1')
    # ]
    best = quads
    res = minimize(optFuncChicane, best, method='nelder-mead', options={'disp': False, 'adaptive': True})
    return res.x

def optFuncVELA(args):
    global dir, lattice
    lattice.modifyElement('EBT-INJ-MAG-QUAD-07', 'k1', args[0])
    lattice.modifyElement('EBT-INJ-MAG-QUAD-08', 'k1', args[1])
    lattice.modifyElement('EBT-INJ-MAG-QUAD-09', 'k1', args[2])
    lattice.modifyElement('EBT-INJ-MAG-QUAD-10', 'k1', args[3])
    lattice.modifyElement('EBT-INJ-MAG-QUAD-11', 'k1', args[4])
    lattice.modifyElement('EBT-INJ-MAG-QUAD-15', 'k1', args[5])
    lattice.track(startfile='VELA', endfile='VELA', track=True)
    twiss.read_sdds_file(dir+'/VELA.twi')
    maxbetaX = max(twiss['betax'])
    maxbetaX = maxbetaX if maxbetaX >= 50 else 0
    maxbetaY = max(twiss['betay'])
    maxbetaY = maxbetaY if maxbetaY >= 50 else 0
    delta = np.sqrt(maxbetaX**2 + maxbetaY**2)
    if delta < 0.1:
        delta = 0.1
    updateOutput('VELA delta = ' + str(delta))
    return delta

def setVELA(quads):
    global dir, lattice
    best = quads
    res = minimize(optFuncVELA, best, method='nelder-mead', options={'disp': False, 'adaptive': True})
    return res.x

def optimise_Lattice():
    global dir, lattice
    dir = './TOMP_SETUP_DBURT'
    lattice = Framework(dir, clean=False, verbose=False)
    lattice.loadSettings('CLA10-BA1_TOMP.def')
    scaling = 5
    if not os.name == 'nt':
        lattice.defineASTRACommand(['mpiexec','-np',str(3*scaling),'/opt/ASTRA/astra_MPICH2.sh'])
        # lattice.defineASTRACommand(['/opt/ASTRA/astra.sh'])
        lattice.defineGeneratorCommand(['/opt/ASTRA/generator.sh'])
        lattice.defineCSRTrackCommand(['/opt/OpenMPI-1.4.3/bin/mpiexec','-n',str(3*scaling),'/opt/CSRTrack/csrtrack_openmpi.sh'])
    lattice.defineElegantCommand(['elegant'])
    # lattice.defineElegantCommand(['mpiexec','-np','6','Pelegant'])
    scaling = 5
    lattice['L01'].file_block['input']['prefix'] = '../basefiles_'+str(scaling)+'/'
    client = unitsC.zmqClient()
    namesK = client.getNamesK(DBURT='CLARA_2_BA1_BA2_2018-11-20-1951.dburt')
    for name, param, k in namesK:
        if param is not None:
            if str(name) in lattice.elementObjects.keys():
                print name, param, k
                lattice.modifyElement(str(name), str(param), float(k))
    # exit()
    # quads = [  21.11058462, -11.36377551,  24.69336696, -22.63264054,  56.07039682, -51.58739658]
    quads = [-13.84327193,  -2.92364597,  18.55449508, -18.53314705]
    quads = setChicane(quads)
    optFuncChicane(quads)
    print 'Chicane = ', quads
    # lattice.modifyElement('EBT-BA1-MAG-QUAD-01', 'k1', -1.22)
    # lattice.modifyElement('EBT-BA1-MAG-QUAD-02', 'k1', 2.18)
    # lattice.modifyElement('EBT-BA1-MAG-QUAD-03', 'k1', -1.22)
    lattice['S02'].file_block['output']['end_element'] = 'EBT-BA1-COFFIN-FOC'

optimise_Lattice()
elegantNoSCOut = open('elegant_NoSC_DBURT_phase_data.csv','wb')
elegantSCOut = open('elegant_SC_DBURT_phase_data.csv','wb')
# ASTRAOut = open('ASTRA_phase_data.csv','wb')
elegentNoSC_csv_out = csv.writer(elegantNoSCOut)
elegentSC_csv_out = csv.writer(elegantSCOut)
# ASTRA_csv_out = csv.writer(ASTRAOut)
for i in range(-25,-24,1):
    set_Phase(i)
    data = track_phase_elegant(i)
    elegentNoSC_csv_out.writerow(data)
    data = track_phase_elegant_SC(i)
    elegentSC_csv_out.writerow(data)
    # data = track_phase_astra(i)
    # ASTRA_csv_out.writerow(data)
