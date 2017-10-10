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
        # print 'SCR_TRACK'
        # First apply a drift through ds/2
        # Next, Calc Y and Y in the middle of the BPM (TP added this)
        self.x = numpy.mean(beam.x)
        self.y = numpy.mean(beam.y)
        self.xSigma = numpy.std(beam.x)
        self.ySigma = numpy.std(beam.y)
        # Finally, apply a second drift through ds/2

        # save
        self.lastTrackedBeam = beam
