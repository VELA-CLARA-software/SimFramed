import read_twiss_file as rtf
import matplotlib.pyplot as plt

twiss = rtf.twiss()
#
twiss.read_astra_emit_files('best_0436'+'/test.in.128.Zemit.128')
print twiss.twiss['alpha_x']
plt.subplot(2, 1, 1)
plt.plot(twiss.twiss['z'], twiss.twiss['beta_x'])
plt.subplot(2, 1, 2)
plt.plot(twiss.twiss['z'], twiss.twiss['cp'])
# plt.plot(beam.slice_bins, 1e8*beam.slice_normalized_horizontal_emittance)
# plt.plot(beam.slice_bins, 1e3*beam.slice_relative_momentum_spread)
plt.grid(True)
plt.show()
