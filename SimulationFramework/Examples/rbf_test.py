import sys, os, math
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname( os.path.abspath(__file__)))))
from SimulationFramework.Framework import *
import SimulationFramework.Modules.read_beam_file as rbf
beam = rbf.beam()

beam.read_HDF5_beam_file('C2V/CLA-C2V-DIA-BPM-01.hdf5', local=True)
print beam.beam['reference_particle']
print beam.beam['starting_position']
