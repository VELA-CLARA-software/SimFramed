# SAMM to Python Conversion
# DJS August 2017
# TP October 2017
# Version 0.2
#
from ComponentBase import ComponentBase
import numpy

class Drift(ComponentBase):
    def __init__(self, length=0, name="", aperture=[]):
        ComponentBase.__init__(self, length, name, aperture)
        self.dummy = 0

    def Track(self, beam):
        # print 'DRIFT_TRACK'
        # Applies the transfer map for a drift to the particles in beam
        # I think the below is NOT as accurate as it might be...
        # this should be updated after some discussions with Wolski
        # The 1 is an approximation for high energy particles...
        d1 = numpy.sqrt(1 - beam.px * beam.px
                        - beam.py * beam.py
                        + 2 * beam.dp / beam.beta
                        + beam.dp * beam.dp)

        # TP: I have kept the 'old' method because I want to compare
        # against SAMM. When We think SAMPL is robust enought we should start
        # making changes

        # abve confirmed!
        # d1 = numpy.sqrt(((beam.gamma*beam.gamma -1)/
        #                   beam.gamma*beam.gamma*beam.beta*beam.beta)
        #                - beam.px * beam.px \
        #                - beam.py * beam.py \
        #                + 2 * beam.dp / beam.beta \
        #                + beam.dp * beam.dp)
        beam.x = beam.x + (self.length * beam.px) / d1
        beam.y = beam.y + self.length * beam.py / d1
        beam.ct = beam.ct + (self.length * (1 - (1 + beam.beta * beam.dp) / d1)
                             / beam.beta)

        # save
        self.lastTrackedBeam = beam
