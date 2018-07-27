import sys, os
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname( os.path.abspath(__file__)))))
import SimulationFramework.Modules.read_beam_file as rbf
import argparse

parser = argparse.ArgumentParser(description='Convert SDDS beam file to ASTRA format')
parser.add_argument('inputfile')
parser.add_argument('outputfile')

args = parser.parse_args()

beam = rbf.beam()
beam.read_pacey_beam_file(args.inputfile)
# beam.beam['longitudinal_reference'] = 't'
beam.write_astra_beam_file(args.outputfile, normaliseZ=4.393091512082e+01)
