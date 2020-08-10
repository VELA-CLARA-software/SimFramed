import sys, os, time
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname( os.path.abspath(__file__)))))
import SimulationFramework.Framework as fw
import SimulationFramework.Modules.read_twiss_file as rtf
import SimulationFramework.Modules.read_beam_file as rbf
import numpy as np

########   Example 1   ########


def example_1():
    # Define a new framework instance, in directory 'C2V'.
    #       "clean" will empty (delete everything!) the directory if true
    #       "verbose" will print a progressbar is true
    lattice = fw.Framework('GPT_Test', clean=False, verbose=False)
    # Load a lattice definition. By default the frameowrk looks in "OnlineModel/MasterLattice/" to find things
    # This example loads a lattice with a CLARA 400Hz gun, and tracking to VELA BA1
    lattice.loadSettings('./CLARA_FE/CLARA_BA1_Gaussian/CLA10-BA1_TOMP_ASTRA.def')
    # This is a scaling parameter
    scaling = 5
    # This defines the number of particles to create at the gun (this is "ASTRA generator" which creates distributions)
    lattice.generator.number_of_particles = 2**(3*scaling)
    # This tracks the beam based on the definitions in the lattice file loaded using "lattice.loadSettings"
    lattice.change_Lattice_Code('All','gpt')
    lattice.track(files=['generator'], preprocess=True, write=True, track=True, postprocess=True)
    for i in np.arange(0.234,0.238,0.0005):
        lattice['CLA-LRG1-GUN-CAV'].phase = -16
        lattice['CLA-LRG1-GUN-SOL'].field_amplitude = i
        lattice.track(files=['injector10'], preprocess=True, write=True, track=True, postprocess=True)
        beam = rbf.beam()
        beam.read_gdf_beam_file('GPT_Test/injector10_out.gdf', position=0.981033)
        print(i, round(1e-6*np.mean(beam.cp),4), round(1e6*beam.normalized_horizontal_emittance,4))
        time.sleep(0.5)
    # print(gdfbeam)
# example_1()
# exit()

########   Example 2   ########
def example_2():
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
    lattice['test_lattice'].write()
    lattice['test_lattice'].run()
    lattice['test_lattice'].postProcess()

# example_2()
# exit()

lattice = fw.Framework('example_ASTRA', clean=False, verbose=True)
lattice.loadSettings('Lattices/clara400_v12_v3.def')
lattice.generator.number_of_particles = 2**(3*4)
lattice.change_Lattice_Code('All','ASTRA', exclude=['injector400','VBC'])
lattice.change_Lattice_Code('VBC','elegant')
lattice.track(files=['generator'])
exit()
# lattice.track(startfile='S04',endfile='S07')

lattice = fw.Framework('example_ASTRA_CSRTrack', clean=False, verbose=True)
lattice.loadSettings('Lattices/clara400_v12_v3.def')
lattice.generator.number_of_particles = 2**(3*4)
lattice.change_Lattice_Code('All','ASTRA', exclude=['injector400','VBC'])
lattice.change_Lattice_Code('VBC','csrtrack')
lattice['S02'].prefix = '../example_ASTRA/'
# lattice.track(startfile='VBC',endfile='S07')

lattice = fw.Framework('example_Elegant', clean=False, verbose=True)
lattice.loadSettings('Lattices/clara400_v12_v3.def')
lattice.generator.number_of_particles = 2**(3*4)
lattice.change_Lattice_Code('All','elegant', exclude=['injector400'])
lattice['S02'].prefix = '../example_ASTRA/'
# lattice.track(startfile='S02', endfile='S07')

lattice = fw.Framework('example_GPT', clean=False, verbose=True)
lattice.loadSettings('Lattices/clara400_v12_v3.def')
lattice.generator.number_of_particles = 2**(3*4)
lattice.change_Lattice_Code('All','elegant', exclude=['injector400'])
lattice.change_Lattice_Code('All','GPT', exclude=['injector400','L02','L03','L04','L4H','VBC','S07'])
lattice['S02'].prefix = '../example_ASTRA/'
# lattice.track(startfile='S07', endfile='S07')
