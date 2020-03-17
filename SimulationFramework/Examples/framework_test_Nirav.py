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
lattice = fw.Framework('Nirav_Test', clean=False, verbose=True)
# Load a lattice definition. By default the frameowrk looks in "OnlineModel/MasterLattice/" to find things
# This example loads a lattice with a CLARA 400Hz gun, and tracking to VELA BA1
lattice.loadSettings('Lattices/clara400_v12_v3_elegant.def')
# This is a scaling parameter npart = 2^(3*scaling)
scaling = 5

# Let's change the VBC angle
print('VBC angle before = ', lattice['bunch_compressor'].angle)
lattice['bunch_compressor'].angle = 0.1185*1.1
print('VBC angle after = ', lattice['bunch_compressor'].angle)
# what about the linacs? Linac1 is in the injector and should be run in ASTRA, so only change if you need to!
print('L02 phase before = ', lattice['CLA-L02-CAV'].phase)
lattice['CLA-L02-CAV'].phase = -18
print('L02 phase after = ', lattice['CLA-L02-CAV'].phase)
# This defines the number of particles to create at the gun (this is "ASTRA generator" which creates distributions)
lattice.generator.number_of_particles = 2**(3*scaling)
# This tracks the beam based on the definitions in the lattice file loaded using "lattice.loadSettings"
# If you have already run the injector before you can start the simulation from, say, the S02
lattice.track(startfile='S02', endfile='S06')
# else run the whole thing
# lattice.track(endfile='S06')
# If you don't want to actually do the tracking you can just write the files out
# lattice.track(startfile='S02', preprocess=False, track=False, postprocess=False)
exit()
