import sys, os
sys.path.append(os.path.abspath(__file__+'/../../../'))
# import SimulationFramework.Framework as fw
# from SimulationFramework.Modules.constraints import *
import SimulationFramework.Modules.read_beam_file as rbf
# import SimulationFramework.Modules.read_twiss_file as rtf
beam = rbf.beam()

etax, etaxp = beam.calculate_gdf_eta('CLARA_FE/CLARA_BA1_Gaussian/GPT/TOMP_GPT_-20/S02_out.gdf')

# beam.read_gdf_beam_file('CLARA_FE/CLARA_BA1_Gaussian/GPT/TOMP_GPT_-20/S02_out.gdf', block=0)
# e, ep = beam.calculate_etax()

print (etax)

s = beam.calculate_gdf_s('CLARA_FE/CLARA_BA1_Gaussian/GPT/TOMP_GPT_-20/S02_emit.gdf')
