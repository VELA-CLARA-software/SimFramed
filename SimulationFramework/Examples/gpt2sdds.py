import sys, os
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname( os.path.abspath(__file__)))))
import SimulationFramework.Modules.read_beam_file as rbf
import argparse

parser = argparse.ArgumentParser(description='Convert GPT GDF beam file to SDDS format')
parser.add_argument('inputfile')
parser.add_argument('outputfile')
parser.add_argument('-r', '--reference', default='z')
parser.add_argument('--info', default=False)
parser.add_argument('-p', '--position', default=None)
parser.add_argument('-t', '--time', default=None)
parser.add_argument('-b', '--block', default=None)
parser.add_argument('-q', '--charge', default=None)
args = parser.parse_args()

beam = rbf.beam()
# beam.read_gdf_beam_file_info(args.inputfile)
if args.info:# or args.info == 'true':
    beam.read_gdf_beam_file_info(args.inputfile)
    exit()
beam.read_gdf_beam_file(args.inputfile, block=args.block, position=args.position, time=args.time, charge=args.charge)
beam.beam['longitudinal_reference'] = args.reference
beam.write_SDDS_file(args.outputfile)
