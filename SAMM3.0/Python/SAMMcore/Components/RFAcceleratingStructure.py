#
from ComponentBase import ComponentBase

from ..SAMMlab import MasterOscillator
from ..SAMMlab import PhysicalConstants as PhyC
import numpy
import math


class RFAcceleratingStructure(ComponentBase):
    def __init__(self, voltage=0, harmonic=1, phase=numpy.pi, length=0,
                 name="", aperture=[], ncell=0, structureType='None'):
        ComponentBase.__init__(self, length, name, aperture)
        # in volts
        self.voltage = voltage
        # frequency (in units of master oscillator frequency)
        self.harmonic = harmonic
        # relative to master oscillator, in radians
        self.phase = phase
        # 1x2 array of elliptical aperture half-axes, in metres
        self.aperture = aperture
        # Number of cells
        self.ncell = ncell
        # Structure type, 'TravellingWave' or 'StandingWave'
        self.structureType = structureType
        # 1 = synchronise with global clock
        self.globalclock = 1

    def setFrequency(self, f):
        MasterOscillator.frequency = f / self.harmonic

    def getFrequency(self):
        return self.harmonic * MasterOscillator.frequency

    def Track(self, beam):
        # beam2 = RFAcceleratingStructure.Track(beam1)
        # Applies the dynamical map for a linac structure to the particles
        # in beam1.  The reference energy is changed by qV*cos(phase).
        mofreq = MasterOscillator.frequency
        nc = self.ncell
        f = self.harmonic * mofreq
        L = self.length / nc
        cosphi = math.cos(self.phase)

        beam.globaltime = (beam.globaltime -
                           math.floor(beam.globaltime * mofreq) / mofreq)

        gt = beam.globaltime * self.globalclock

        dE = beam.species.charge * self.voltage * cosphi / nc
        E0 = beam.energy
        E1 = E0 + dE

        if abs(dE / E0) > 1e-6:
            if self.structuretype == 'TravellingWave':
                logE = math.log(1 + dE / E0)
                r11 = 1 - logE / 2
                r12 = L * logE * E0 / dE
                r21 = -dE * logE / 4 / L / E1
                r22 = E0 * (2 + logE) / 2 / E1
            elif self.structuretype == 'StandingWave':
                L1 = PhyC.SpeedOfLight / 2 / f
                if abs(L / L1 - 1) > 1e-2:
                    print('RFAcceleratingStructure:BadLength  ',
                          'RFAcceleratingStructure.length should be c/2f.')
            else:
                print('RFAcceleratingStructure:UnrecognisedType  ',
                      'RFAcceleratingStructure.structuretype should',
                      ' be StandingWave or TravellingWave.')
        else:
            r11 = 1
            r12 = L * (1 - dE / 2 / E0)
            r21 = 0
            r22 = 1 - dE / E0
        x0 = beam.x
        px0 = beam.px
        y0 = beam.y
        py0 = beam.py
        ct0 = beam.ct
        dp0 = beam.dp
        for n in range(1, nc + 1):
            print 'hi'
            # First, apply a drift map through L/2
            # to the longitudinal coordinate
            beta0 = beam.__beta
            d1 = math.sqrt(1 - px0 * px0 - py0 * py0 +
                           2 * dp0 / beta0 + dp0 * dp0)
            ct0 = ct0 + L * (1 - (1 + beta0 * dp0) / d1) / beta0 / 2

            # Now apply the RF structure map to the transverse variables
            # and the momentum deviation
            x1 = r11 * x0 + r12 * px0
            px0 = r21 * x0 + r22 * px0
            x0 = x1

            y1 = r11 * y0 + r12 * py0
            py0 = r21 * y0 + r22 * py0
            y0 = y1

            P0 = beam.momentum
            Edc = (dp0 + 1 / beta0) * P0

            beam.energy = E1
            beta0 = beam.beta
            P0 = beam.momentum

            t = gt - ct0 / PhyC.SpeedOfLight
            Edc = Edc +
            beam.species.charge * self.voltage *
            math.cos(2 * math.pi * f * t + self.phase) / nc / PhyC.SpeedOfLight

            dp0 = Edc / P0 - 1 / beta0

            # Finally, apply a drift map through L/2
            # to the longitudinal coordinate
            d1 = math.sqrt(1 - px0 * px0 - py0 * py0 +
                           2 * dp0 / beta0 + dp0 * dp0)
            ct0 = ct0 + L * (1 - (1 + beta0 * dp0) / d1) / beta0 / 2
        beam.x = x0
        beam.px = px0
        beam.y = y0
        beam.py = py0
        beam.ct = ct0
        beam.dp = dp0

        # save
        self.lastTrackedBeam = beam
