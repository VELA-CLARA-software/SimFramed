from sourceCode.SAMPLcore.SAMPLlab import Beam
from sourceCode.SAMPLcore.Particles import Electron
from sourceCode.SAMPLcore.Particles import Positron
# from SAMMcore.Particles import Proton
from sourceCode.SAMPLcore.SAMPLlab import PhysicalUnits
import numpy as np


class createBeam():

    def __init__(self, PIL_Ctrl=None):
        self.PIL_Ctrl = PIL_Ctrl

    def guassian(self, x=0.0, y=0.0, sigmaX=0.001, sigmaY=0.001,
                 particle='Electron', number=1000,
                 Energy=4.5 * PhysicalUnits.MeV, charge=250e-9):
        if particle == 'Electron':
            beam = Beam.Beam(species=Electron.Electron, energy=Energy)
        if particle == 'Positron':
            beam = Beam.Beam(species=Positron.Positron, energy=Energy,
                             bunchcharge=charge)
        # if particle == 'Proton':
        #     beam = Beam.Beam(species=Proton.Proton, energy = Energy)
        # Generate particles
        xArray = np.random.normal(loc=x, scale=sigmaX, size=number)
        yArray = np.random.normal(loc=y, scale=sigmaY, size=number)
        ctArray = np.zeros(number)
        pxArray = np.zeros(number)
        pyArray = np.zeros(number)
        dpArray = np.zeros(number)

        p = np.array([xArray, pxArray, yArray, pyArray, ctArray, dpArray])
        p = p.transpose()
        beam.particles = p

        return beam
