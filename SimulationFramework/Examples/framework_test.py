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

lattice = fw.Framework('example', clean=False, verbose=True)
lattice.loadSettings('FEBE/FEBE.def')
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
lattice['S02'].prefix = '../CLARA/basefiles_5/'
lattice.change_Lattice_Code('All','elegant')
# lattice.track(startfile='S02')
lattice.read_Lattice('test_lattice', {'code': 'elegant', 'charge': {'cathode': False, 'space_charge_mode': '3D'}, 'input': {}, 'output': {'start_element': 'CLA-S07-MARK-01', 'end_element': 'CLA-FEB-W-START-01'}})
lattice['test_lattice'].preProcess()
lattice['test_lattice'].optimisation = fw.elegantOptimisation(lattice['test_lattice'], variables={'test': {'item': 'betax', 'lower': 1, 'upper': 10}})
print lattice['test_lattice'].optimisation.commandObjects['test'].write()
exit()
lattice['test_lattice'].write()
lattice['test_lattice'].run()
lattice['test_lattice'].postProcess()
exit()


lattice['S02'].prefix = '../CLARA/basefiles_5/'
# lattice.setSubDirectory('example_astra')
# lattice.track(startfile='FMS', endfile='FMS')
lattice.setSubDirectory('example_elegant')
lattice.change_Lattice_Code('All','elegant')
lattice.load_changes_file('./CLARA/Elegant_Genesis/Short-240/transverse_best_changes.yaml')
# lattice.modifyElement('CLA-S07-DCP-01', 'scale_kick', 1e9)
lattice.track(startfile='S02', endfile='FMS')
# # lattice.track(startfile='S02')#,preprocess=True, track=False, postprocess=False)
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
