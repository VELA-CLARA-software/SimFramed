import read_beam_file as raf
import matplotlib.pyplot as plt

beam = raf.beam()
#
beam.slice_length = 0.025e-12
#
''' Testing reading ASTRA file '''
beam.read_astra_beam_file("F:\My Documents\Mathematica\CLARA\CLARA_V12\Short_240_10k_noTDC\clara.4636.018")
print beam.cpz/(2.998*10**8)*1.6*10**-19
''' Testing reading GDF file with pages '''
# gdfbeam = beam.read_gdf_beam_file_object('gdf_test.gdf')
# beam.read_gdf_beam_file(gdfbeam, block=0)
# positions = gdfbeam.positions
# beam.read_gdf_beam_file(gdfbeam, position=positions[1])

''' Testing reading Elegant SDDS file with pages '''
# beam.read_sdds_file('W-END.sdds', charge=250e-12)
# plt.figure()
# plt.plot(beam.slice_bins, abs(beam.slice_peak_current))
# plt.plot(beam.slice_bins, 1e8*beam.slice_normalized_horizontal_emittance)
# plt.plot(beam.slice_bins, 1e3*beam.slice_relative_momentum_spread)
# plt.grid(True)
# plt.show()
# print beam.twiss_analysis
# print beam.twiss_analysis_corrected
# print beam.horizontal_emittance
# print beam.horizontal_emittance_corrected
# beam.write_astra_beam_file('scapaS1.0150.001', normaliseZ=0)
beam.write_vsim_beam_file('clara.4636.018.vsim', normaliseT=True)
