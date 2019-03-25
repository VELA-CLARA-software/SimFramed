import sys, os
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname( os.path.abspath(__file__)))))
import SimulationFramework.Framework as fw
import SimulationFramework.Modules.read_twiss_file as rtf
import numpy as np

# lattice = framework('C2V', clean=False, verbose=False)
# lattice.loadSettings('Lattices/cla400-ba1.def')
# if not os.name == 'nt':
#     scaling = 5
#     lattice.defineASTRACommand(['mpiexec','-np',str(3*scaling),'/opt/ASTRA/astra_MPICH2.sh'])
#     lattice.defineGeneratorCommand(['/opt/ASTRA/generator.sh'])
#     lattice.defineCSRTrackCommand(['/opt/OpenMPI-1.4.3/bin/mpiexec','-n',str(3*scaling),'/opt/CSRTrack/csrtrack_openmpi.sh'])
#     lattice.generator.number_of_particles = 2**(3*scaling)
# else:
#     lattice.generator.number_of_particles = 2**(3*3)
# lattice.track()#startfile='C2V')

lattice = fw.Framework('example', clean=True, verbose=True)
lattice.loadSettings('Lattices/clara400_v12_v3.def')
if not os.name == 'nt':
    scaling = 5
    lattice.defineASTRACommand(['mpiexec','-np',str(3*scaling),'/opt/ASTRA/astra_MPICH2.sh'])
    # lattice.defineASTRACommand(['/opt/ASTRA/astra.sh'])
    lattice.defineGeneratorCommand(['/opt/ASTRA/generator.sh'])
    lattice.defineCSRTrackCommand(['/opt/OpenMPI-1.4.3/bin/mpiexec','-n',str(3*scaling),'/opt/CSRTrack/csrtrack_openmpi.sh'])
    lattice.generator.number_of_particles = 2**(3*scaling)
else:
    lattice.generator.number_of_particles = 2**(3*5)
lattice.defineElegantCommand(['elegant'])


lattice['S02'].prefix = '../CLARA/basefiles_4/'
oldl02grad = lattice.getElement('CLA-L02-CAV', 'field_amplitude')
oldl02phase = lattice.getElement('CLA-L02-CAV', 'phase')
oldl03grad = lattice.getElement('CLA-L03-CAV', 'field_amplitude')
oldl03phase = lattice.getElement('CLA-L03-CAV', 'phase')
oldl04grad = lattice.getElement('CLA-L04-CAV', 'field_amplitude')
oldl04phase = lattice.getElement('CLA-L04-CAV', 'phase')
oldl0xgrad = lattice.getElement('CLA-L4H-CAV', 'field_amplitude')
oldl0xphase = lattice.getElement('CLA-L4H-CAV', 'phase')
lattice.modifyElement('CLA-L02-CAV', 'field_amplitude', np.sqrt(2)*11.25e6 )
lattice.modifyElement('CLA-L03-CAV', 'field_amplitude', np.sqrt(2)*11.25e6 )
lattice.modifyElement('CLA-L04-CAV', 'field_amplitude', np.sqrt(2)*12e6 )
lattice.modifyElement('CLA-L4H-CAV', 'field_amplitude', np.sqrt(2)*15.75e6 )
lattice.modifyElement('CLA-L02-CAV', 'phase', oldl02phase + 5.0)
fw.save_change_file(lattice)
print lattice['S02'].prefix
print lattice.S02['CLA-S02-MAG-QUAD-01']['k1l']
lattice.modifyElement('CLA-S02-MAG-QUAD-01','k1', 10)
print lattice['S02'].getElement('CLA-S02-MAG-QUAD-01','k1')
lattice['CLA-S02-MAG-QUAD-01'].k1 = 13.15
print lattice['S02'].getElement('CLA-S02-MAG-QUAD-01','k1')
# lattice.track(files=['generator', 'injector400'], track=True, postprocess=True)
lattice.change_Lattice_Code('All','elegant')
# lattice.track(startfile='S02')#,preprocess=True, track=False, postprocess=False)
# lattice['bunch_compressor'].set_angle(-0.0)
# lattice.track(startfile='S02', endfile='S07')#,preprocess=True, track=False, postprocess=False)

# exit()
# def f(ph):
#     lattice.setSubDirectory('example/'+'L01_phase_'+str(ph))
#     lattice['CLA-L01-CAV'].phase = ph
#     lattice.track()
#
# from multiprocessing import Pool
# if __name__ == '__main__':
#     p = Pool(10)
#     startphase = lattice['CLA-L01-CAV'].phase
#     vals = startphase + np.arange(-10,10.01,5)
#     # print vals
#     print(p.map(f, vals))
