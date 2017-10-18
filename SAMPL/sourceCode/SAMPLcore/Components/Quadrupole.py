# SAM to Python Conversion
# DJS August 2017
# Version 0.1
#
from ComponentBase import ComponentBase
import numpy
from ..SAMPLlab import Utilities
import Drift


class Quadrupole(ComponentBase):
    def __init__(self, length=0, name="", aperture=[], gradient=0):
        ComponentBase.__init__(self, length, name, aperture)

        # quadrupole gradient, in tesla/metre
        self.gradient = gradient
        # self.__bx = None
        # self.__by = None
        # self.__bz = None
        self.k1 = None
        self.absk1 = None
        self.ds = self.length

    def Track(self, beam):
        # normalised gradient
        self.k1 = numpy.divide(self.gradient, beam.rigidity)
        self.absk1 = abs(self.k1)
        # print('Quad = ',self.gradient,beam.rigidity)
        if self.k1 > 0:
            self.TrackFQuad(beam)
        elif self.k1 < 0:
            self.TrackDQuad(beam)
        else:
            # print('Quad TrackDrift')
            self.TrackDrift(beam)

        # save
        self.lastTrackedBeam = beam

    def TrackFQuad(self, beam):
        # print 'QUAD_TRACK'
        d1 = numpy.sqrt(1 + numpy.divide(2 * beam.dp, beam.beta) +
                        beam.dp * beam.dp)
        w = numpy.sqrt(numpy.divide(self.k1, d1))

        xs = numpy.sin(w * self.ds)
        xc = numpy.cos(w * self.ds)
        ys = numpy.sinh(w * self.ds)
        yc = numpy.cosh(w * self.ds)
        xs2 = numpy.sin(2 * w * self.ds)
        ys2 = numpy.sinh(2 * w * self.ds)

        d0 = 1 / beam.beta + beam.dp
        d2 = -d0 / d1 / d1 / d1 / 2

        c0 = (1 / beam.beta - d0 / d1) * self.ds
        c11 = self.k1 * self.k1 * d2 * (xs2 / w - 2 * self.ds) / w / w / 4
        c12 = -self.k1 * d2 * xs * xs / w / w
        c22 = d2 * (xs2 / w + 2 * self.ds) / 4
        c33 = self.k1 * self.k1 * d2 * (ys2 / w - 2 * self.ds) / w / w / 4
        c34 = self.k1 * d2 * ys * ys / w / w
        c44 = d2 * (ys2 / w + 2 * self.ds) / 4
        beam.ct = (beam.ct + c0
                   + c11 * beam.x * beam.x
                   + c12 * beam.x * beam.px
                   + c22 * beam.px * beam.px
                   + c33 * beam.y * beam.y
                   + c34 * beam.y * beam.py
                   + c44 * beam.py * beam.py)

        x0 = beam.x
        beam.x = beam.x * xc + beam.px * (xs * w / self.absk1)
        beam.px = -x0 * self.k1 * xs / w + beam.px * xc
        y0 = beam.y
        beam.y = beam.y * yc + beam.py * (ys * w / self.absk1)
        beam.py = y0 * self.k1 * ys / w + beam.py * yc

    def TrackDQuad(self, beam):
        # print 'QUAD_TRACK'
        d1 = numpy.sqrt(1 + 2 * numpy.divide(beam.dp, beam.beta) + numpy.multiply(beam.dp, beam.dp))
        w = numpy.sqrt(self.absk1 / d1)

        xs = numpy.sinh(w * self.ds)
        xc = numpy.cosh(w * self.ds)
        ys = numpy.sin(w * self.ds)
        yc = numpy.cos(w * self.ds)
        xs2 = numpy.sinh(2 * w * self.ds)
        ys2 = numpy.sin(2 * w * self.ds)

        d0 = 1 / beam.beta + beam.dp
        d2 = -d0 / d1 / d1 / d1 / 2
        c0 = (1 / beam.beta - d0 / d1) * self.ds
        c11 = self.k1 * self.k1 * d2 * (xs2 / w - 2 * self.ds) / w / w / 4
        c12 = -self.k1 * d2 * xs * xs / w / w
        c22 = d2 * (xs2 / w + 2 * self.ds) / 4
        c33 = self.k1 * self.k1 * d2 * (ys2 / w - 2 * self.ds) / w / w / 4
        c34 = self.k1 * d2 * ys * ys / w / w
        c44 = d2 * (ys2 / w + 2 * self.ds) / 4
        beam.ct = (beam.ct + c0
                   + c11 * beam.x * beam.x
                   + c12 * beam.x * beam.px
                   + c22 * beam.px * beam.px
                   + c33 * beam.y * beam.y
                   + c34 * beam.y * beam.py
                   + c44 * beam.py * beam.py)
        x0 = beam.x
        beam.x = beam.x * xc + beam.px * (xs * w / self.absk1)
        beam.px = -x0 * self.k1 * xs / w + beam.px * xc
        y0 = beam.y
        beam.y = beam.y * yc + beam.py * (ys * w / self.absk1)
        beam.py = y0 * self.k1 * ys / w + beam.py * yc

    def TrackDrift(self, beam):
        # print 'QUAD_TRACK'
        drift = Drift.Drift(length=self.length)
        drift.Track(beam)

    def GetBField(self, beam):
        # [bx, by, bz] = Quadrupole.GetBField(beam)
        # Returns the magnetic field (in tesla) at the locations of
        # the particles in the beam.
        bx = beam.y * self.gradient
        by = beam.x * self.gradient
        bz = [0.0] * len(beam.y)
        return [bx, by, bz]


    def TrackSpin(self, beam):
        # Quadrupole.Trackspin(beam)
        # Tracks particle spins through a Quadrupole.
            [bx, by, bz] = self.GetBField(beam)
            Utilities.SpinRotation(beam,bx,by,bz,self.length)
