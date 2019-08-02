import sys, os
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname( os.path.abspath(__file__)))))
import SimulationFramework.Framework as fw
import SimulationFramework.Modules.read_twiss_file as rtf
import SimulationFramework.Modules.read_beam_file as rbf
import numpy as np

########   Example 1   ########

# Define a new framework instance, in directory 'C2V'.
#       "clean" will empty (delete everything!) the directory if true
#       "verbose" will print a progressbar is true
lattice = fw.Framework('GPT_Test', clean=False, verbose=False)
# Load a lattice definition. By default the frameowrk looks in "OnlineModel/MasterLattice/" to find things
# This example loads a lattice with a CLARA 400Hz gun, and tracking to VELA BA1
lattice.loadSettings('./CLARA_FE/CLARA_BA1_Gaussian/CLA10-BA1_TOMP_ASTRA.def')
# This is a scaling parameter 
scaling = 5
# This defines the location of tracking codes - On windows this uses versions in "OnlineModel/MasterLattice/Codes",
# but you will need to define these yourself on linux
lattice.defineASTRACommand(location=['/opt/ASTRA/astra_MPICH2.sh'])
lattice.defineGeneratorCommand(location=['/opt/ASTRA/generator.sh'])
lattice.defineCSRTrackCommand(location=['/opt/CSRTrack/csrtrack.sh'])
# This defines the number of particles to create at the gun (this is "ASTRA generator" which creates distributions)
lattice.generator.number_of_particles = 2**(3*scaling)
# This tracks the beam based on the definitions in the lattice file loaded using "lattice.loadSettings"
lattice.change_Lattice_Code('All','gpt')
lattice['S02'].prefix = '../CLARA_FE/CLARA_BA1_Gaussian/TOMP_SETUP_-10/'
lattice['S02'].sample_interval = 2**(3*1)
lattice.track(files=['S02', 'C2V'], preprocess=True, write=True, track=True, postprocess=True)
# beam = rbf.beam()
# beam.read_gdf_beam_file('GPT_Test/S02_out.gdf', position=5.72687)
# print(beam.zn)
# print(gdfbeam)
exit()

########   Example 2   ########
lattice = fw.Framework('example', clean=False, verbose=True)
lattice.loadSettings('FEBE/FEBE.def')
if not os.name == 'nt': # i.e. we are on linux
    scaling = 5
    lattice.defineASTRACommand(scaling=scaling)
    # lattice.defineASTRACommand(['/opt/ASTRA/astra.sh'])
    lattice.defineCSRTrackCommand(scaling=scaling)
    lattice.generator.number_of_particles = 2**(3*scaling)
else:
    lattice.generator.number_of_particles = 2**(3*5)

lattice['S02'].prefix = '../CLARA/basefiles_5/'
lattice.change_Lattice_Code('All','elegant')
# lattice.track(startfile='S02')
lattice.read_Lattice('test_lattice', {'code': 'elegant', 'charge': {'cathode': False, 'space_charge_mode': '3D'}, 'input': {}, 'output': {'start_element': 'CLA-S07-MARK-01', 'end_element': 'CLA-FEB-W-START-01'}})
lattice['test_lattice'].preProcess()
lattice['test_lattice'].optimisation = fw.elegantOptimisation(lattice['test_lattice'], variables={'test': {'item': 'betax', 'lower': 1, 'upper': 10}})
print (lattice['test_lattice'].optimisation.commandObjects['test'].write())
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
