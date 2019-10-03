import sys, os
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/../../../../')
import SimulationFramework.Framework as fw
import SimulationFramework.Modules.read_twiss_file as rtf
import SimulationFramework.Modules.read_beam_file as rbf

beam = rbf.beam()
beam.read_gdf_beam_file('GPT_CSR/TOMP_GPT_CSR_-20/S02_out.gdf', position=0.5210126725216484)
