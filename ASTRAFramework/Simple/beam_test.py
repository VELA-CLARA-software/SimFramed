import read_beam_file as raf
import matplotlib.pyplot as plt

beam = raf.beam()
#
beam.slice_length = 0.025e-12
#
beam.read_astra_beam_file('best_0436'+'/test.in.128.4929.128')
# plt.figure()
# plt.plot(beam.slice_bins, abs(beam.slice_peak_current))
# plt.plot(beam.slice_bins, 1e8*beam.slice_normalized_horizontal_emittance)
# plt.plot(beam.slice_bins, 1e3*beam.slice_relative_momentum_spread)
# plt.grid(True)
# plt.show()
print beam.sliceAnalysis()
# #
beam.read_sdds_file('best_0436'+'/test.in.128.4929.128.sdds', charge=250e-12)
print beam.sliceAnalysis()
#
beam.read_gdf_beam_file('best_0436'+'/test.in.128.4929.128.gdf', charge=250e-12)
print beam.sliceAnalysis()
#
print beam.twiss_analysis
# print beam.twiss_analysis_corrected
print beam.horizontal_emittance_90
