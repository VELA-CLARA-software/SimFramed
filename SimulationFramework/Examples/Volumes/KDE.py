import sys, os, time
sys.path.append(os.path.abspath(__file__+'/../../../../'))
import SimulationFramework.Modules.read_beam_file as rbf
from fastkde import fastKDE
import numpy as np

def measure(n):
    "Measurement model, return two coupled measurements."
    m1 = np.random.normal(size=n)
    m2 = np.random.normal(scale=0.5, size=n)
    return m1+m2, m1-m2

# m1, m2 = measure(2000)
start = time.clock()
dirname = '../CLARA/CLARA_best_longitudinal'
beam = rbf.beam()
beam.read_HDF5_beam_file(dirname+'/CLA-FMS-APER-01.hdf5')
print 'Read beam files = ', time.clock() - start
m1 = 100*beam.x
m2 = 100*beam.y
print max(beam.x)
xmin = -0.03
xmax = 0.03
ymin = -0.03
ymax = 0.03

X, Y = np.mgrid[xmin:xmax:100j, ymin:ymax:100j]
print 'Grid = ', time.clock() - start
positions = np.vstack([X.ravel(), Y.ravel()])
print 'Positions = ', time.clock() - start
values = [m1, m2]
print 'Values = ', time.clock() - start
kernel,axes = fastKDE.pdf(m1, m2)
print 'Kernel = ', time.clock() - start
kpos = kernel(positions)
print 'Kpos = ', time.clock() - start
Z = np.reshape(kpos.T, X.shape)
print 'Z= ', time.clock() - start
import matplotlib.pyplot as plt
fig, ax = plt.subplots()
ax.imshow(np.rot90(Z), cmap=plt.cm.gist_earth_r,
          extent=[xmin, xmax, ymin, ymax])
# ax.plot(m1[::10], m2[::10], 'k.', markersize=2)
ax.set_xlim([xmin, xmax])
ax.set_ylim([ymin, ymax])
plt.show()
