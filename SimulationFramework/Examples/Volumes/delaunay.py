import sys, os, time, math
sys.path.append(os.path.abspath(__file__+'/../../../../'))
import SimulationFramework.Modules.read_beam_file as rbf
from scipy import stats
import numpy as np
from scipy.spatial import ConvexHull, Delaunay, Voronoi, voronoi_plot_2d
import matplotlib.pyplot as plt

start = time.clock()
dirname = '../CLARA/CLARA_best_longitudinal'
beam = rbf.beam()
beam.read_HDF5_beam_file(dirname+'/CLA-FMS-APER-01.hdf5')
print 'Read beam files = ', time.clock() - start

x_test = beam.x[::4096]
y_test = beam.y[::4096]
points = np.array([x_test, y_test]).T

def PolyArea(tri):
    x,y = tri.T
    return 0.5*np.abs(np.dot(x,np.roll(y,1))-np.dot(y,np.roll(x,1)))

vor = Voronoi(points)
print vor.regions
voronoi_plot_2d(vor)
# regions = [vor.vertices[vor.regions[x]] for x in vor.point_region]
# z_test = np.array([ConvexHull(x).volume for x in regions])
# z_test = z_test / max(z_test)
# print max(z_test)
# import matplotlib.pyplot as plt
#
# plt.scatter(x_test, y_test, c=z_test, s=1)
#
plt.show()
