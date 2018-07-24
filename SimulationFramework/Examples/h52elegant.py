import sys, os
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname( os.path.abspath(__file__)))))
import SimulationFramework.Modules.read_beam_file as rbf
import argparse

parser = argparse.ArgumentParser(description='Convert hdf5 beam file to SDDS format')
parser.add_argument('inputfile')
parser.add_argument('outputfile')

args = parser.parse_args()

beam = rbf.beam()
beam.read_HDF5_beam_file(args.inputfile)
beam.write_SDDS_file(args.outputfile)
