import sys, os
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname( os.path.abspath(__file__)))))
import SimulationFramework.Modules.read_beam_file as rbf
import argparse

parser = argparse.ArgumentParser(description='Convert GPT GDF beam file to ASTRA format')
parser.add_argument('inputfile')
parser.add_argument('outputfile')

args = parser.parse_args()

beam = rbf.beam()
# beam.read_gdf_beam_file_info(args.inputfile)
beam.read_gdf_beam_file(args.inputfile, block=0)
beam.beam['longitudinal_reference'] = 'z'
beam.write_astra_beam_file(args.outputfile, normaliseZ=True)
