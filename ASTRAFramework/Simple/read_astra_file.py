import numpy as np
from  scipy.constants import *

def normalise_to_ref_particle(array, index=0,subtractmean=False):
    array[1:] = array[0] + array[1:]
    if subtractmean:
        array = array - np.mean(array)
    return array

def read_astra_beam_file(file):
    beam = {}
    x, y, z, px, py, pz, clock, charge, index, status = np.loadtxt(file, unpack=True)
    z = normalise_to_ref_particle(z, subtractmean=True)
    t = -1*z / speed_of_light
    pz = normalise_to_ref_particle(pz, subtractmean=False)
    p = np.sqrt(px**2 + py**2 + pz**2)
    #print p
    beam['x'] = x
    beam['y'] = y
    beam['z'] = z
    beam['t'] = t
    beam['px'] = px
    beam['py'] = py
    beam['pz'] = pz
    beam['p'] = p
    beam['xp'] = px/pz
    beam['yp'] = py/pz
    beam['clock'] = clock
    beam['charge'] = charge
    beam['index'] = index
    beam['status'] = status
    return beam
