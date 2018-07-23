import sys, os
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname( os.path.abspath(__file__)))))
import SimulationFramework.Modules.read_beam_file as rbf
import argparse

parser = argparse.ArgumentParser(description='Rematch Twiss parameters for an astra beam file')
parser.add_argument('inputfile')
parser.add_argument('outputfile')
parser.add_argument('-betax', default=False, type=float)
parser.add_argument('-alphax', default=False, type=float)
parser.add_argument('-betay', default=False, type=float)
parser.add_argument('-alphay', default=False, type=float)
parser.add_argument('-sigmax', default=False, type=float)
parser.add_argument('-sigmay', default=False, type=float)
args = parser.parse_args()

beam = rbf.beam()
beam.read_astra_beam_file(args.inputfile)
if args.sigmax is not False:
    args.betax = (beam.beta_x/((beam.Sx)/args.sigmax)**2)
if args.sigmay is not False:
    args.betay = (beam.beta_y/((beam.Sy)/args.sigmay)**2)
beam.rematchXPlane(beta=args.betax, alpha=args.alphax)
beam.rematchYPlane(beta=args.betay, alpha=args.alphay)
beam.write_astra_beam_file(args.outputfile)
beam.write_SDDS_file(args.outputfile+'.sdds')
