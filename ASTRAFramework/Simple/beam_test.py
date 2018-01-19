import read_beam_file as raf
import matplotlib.pyplot as plt

beam = raf.beam()
#
beam.slice_length = 0.025e-12
#
''' Testing reading ASTRA file '''
beam.read_astra_beam_file('best_0955/test.in.128.4929.128')
import numpy as np
sigmat = 1e12*np.std(beam.t)
sigmap = np.std(beam.p)
meanp = np.mean(beam.p)
emitx = 1e6*beam.normalized_horizontal_emittance
emity = 1e6*beam.normalized_vertical_emittance
fitp = 100*sigmap/meanp
fhcfield = 58
peakI, peakIMomentumSpread, peakIEmittanceX, peakIEmittanceY, peakIMomentum = beam.sliceAnalysis()
constraintsList = {
    'peakI': {'type': 'greaterthan', 'value': abs(peakI), 'limit': 400, 'weight': 20},
    'peakIMomentumSpread': {'type': 'lessthan', 'value': peakIMomentumSpread, 'limit': 0.1, 'weight': 3},
    'peakIEmittanceX': {'type': 'lessthan', 'value': 1e6*peakIEmittanceX, 'limit': 0.25, 'weight': 2.5},
    'peakIEmittanceY': {'type': 'lessthan', 'value': 1e6*peakIEmittanceY, 'limit': 0.25, 'weight': 2.5},
    'peakIMomentum': {'type': 'greaterthan','value': 1e-6*peakIMomentum, 'limit': 210, 'weight': 10},
    'linac fields': {'type': 'lessthan', 'value': max([1]), 'limit': 32, 'weight': 100},
    '4hc field': {'type': 'lessthan', 'value': fhcfield, 'limit': 35, 'weight': 100},
    'horizontal emittance': {'type': 'lessthan', 'value': emitx, 'limit': 1, 'weight': 3},
    'vertical emittance': {'type': 'lessthan', 'value': emity, 'limit': 1, 'weight': 3},
}
print constraintsList
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
# beam.write_vsim_beam_file('W-END.vsim', normaliseT=True)
