# SAM to Python Conversion
# DJS August 2017
# Version 0.1
#
import numpy
import math
from ComponentBase import ComponentBase


class Screen(ComponentBase):
    def __init__(self, name="", aperture=[]):
        # Could use the apertures here for the slits that we
        # have on screen conveyer belts.
        # The length will ALWAYS be 0.0m
        length = 0.0
        ComponentBase.__init__(self, length, name, aperture)
        # super(ComponentBase, self).__init__(length, name, aperture)
        self.dummy = 0
        self.x = None
        self.y = None
        self.xSigma = None
        self.ySigma = None
    def Track(self, beam):
        # First apply a drift through ds/2
        d1 = math.sqrt(1 - beam.px**2 - beam.py**2 + 2 * beam.dp / beta0 +
                       beam.dp**2)
        x1 = beam.x + ds * beam.px / d1 / 2
        y1 = beam.y + ds * beam.py / d1 / 2
        ct1 = beam.ct + ds * (1 - (1 + beam.dp * beta0) / d1) / beta0 / 2
        # Next, Calc Y and Y in the middle of the BPM (TP added this)
        self.x = numpy.mean(x1)
        self.y = numpy.mean(y1)
        self.xSigma = numpy.std(x1)
        self.ySigma = numpy.std(y1)
        # Finally, apply a second drift through ds/2
