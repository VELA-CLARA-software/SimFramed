import sys, os
import numpy as np
sys.path.append(os.path.abspath(__file__+'/../../../../'))
import SimulationFramework.Modules.read_beam_file as rbf
from scipy import spatial
from scipy.special import gamma
import timeit
from functools import partial


def kNN_Volume(data):
    k = int(np.floor(np.sqrt(len(data))))
    tree = spatial.cKDTree(data, balanced_tree=False)
    d = np.shape(data)[-1]

    distance, i = tree.query(data, k, n_jobs=-1)
    radii = distance[:,-1]**d
    dens = (k * gamma(d/2 + 1)) / (np.pi**(d/2)) / radii

    xy = data[np.argsort(dens)][-int(np.round(0.2*len(data))):]
    return spatial.ConvexHull(xy, qhull_options='QJ').volume

def kNN_Volume_Beam(beam):
    data = np.array(zip(beam.x, beam.y, beam.z, beam.cpx/beam.cp, beam.cpy/beam.cp, beam.cpz/beam.cp))
    return kNN_Volume(data)

def kNN_Volume_File(file):
    beam = rbf.beam()
    beam.read_HDF5_beam_file(file)
    return kNN_Volume_Beam(beam)

if __name__ == "__main__":

    dirname = '../CLARA/CLARA_best_longitudinal'
    print kNN_Volume_File(dirname+'/CLA-FMS-APER-01.hdf5')
