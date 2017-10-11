# SAM to Python Conversion
# DJS August 2017
# Version 0.1
#
from ComponentBase import ComponentBase
from ..SAMMlab import Beam
import numpy

class Drift(ComponentBase):
    def __init__(self, length=0, name="", aperture=[]):
        ComponentBase.__init__(self, length, name, aperture)
        # super(ComponentBase, self).__init__(length, name, aperture)
        self.dummy = 0

    def Track(self, beam):
        # print 'DRIFT_TRACK'
        # print 'Tracking - DRIFT:  particles'
        # Applies the transfer map for a drift to the particles in beam
        #
        # I think the below is NOT as accurate as it might be...
        # this should be updated after some discussions with Wolski
        # The 1 is an approximation for high energy particles...
        d1 = numpy.sqrt(1 - beam.px * beam.px
                        - beam.py * beam.py
                        + 2 * beam.dp / beam.beta
                        + beam.dp * beam.dp)
        # abve confirmed!
        # d1 = numpy.sqrt((beam.gamma*beam.gamma -1)/beam.gamma*beam.gamma*beam.beta*beam.beta
        #                - beam.px * beam.px \
        #                - beam.py * beam.py \
        #                + 2 * beam.dp / beam.beta \
        #                + beam.dp * beam.dp)

        # print('numpy.divide( beam.px, d1) = ', numpy.divide(beam.px, d1))
        # remember these are relative to the reference particle (!)
        # beam.x  = beam.x  + self.length*numpy.divide(beam.px,d1)
        beam.x = beam.x + (self.length * beam.px) / d1
        # beam.y  = beam.y  + self.length*numpy.divide(beam.py,d1)
        beam.y = beam.y + self.length * beam.py / d1

        # beam.ct = beam.ct + self.length*numpy.divide\
        #         ((1 - numpy.divide(1 + beam.beta*beam.dp,d1)),beam.beta)

        beam.ct = beam.ct + self.length * (1 - (1 + beam.beta * beam.dp) / d1) / beam.beta

        # save
        self.lastTrackedBeam = beam



# classdef Drift < handle
#     % Drift class
#     %
#     % Properties:
#     %   name
#     %   length
#     %   aperture
#     %
#     % Methods:
#     %   Track
#
#     properties
#         name     = ''; % string
#         length   = 0;  % in metres
#         aperture = []; % 1x2 array of elliptical aperture half-axes, in metres
#     end % properties
#
#     methods
#
#         function beam = Track(drift,beam)
#             % beam2 = Drift.Track(beam1)
#             % Applies the transfer map for a drift to the particles
#             % in beam1.
#
#             [x0, px0, y0, py0, ct0, dp0] = beam.GetParticles();
#
#             beta0  = beam.beta;
#
#             ds  = drift.length;
#
#             d1  = sqrt(1 - px0.*px0 - py0.*py0 + 2*dp0/beta0 + dp0.*dp0);
#
#             x1  = x0  + ds*px0./d1;
#             y1  = y0  + ds*py0./d1;
#             ct1 = ct0 + ds*(1 - (1 + beta0*dp0)./d1)/beta0;
#
#             beam.SetParticles(x1, px0, y1, py0, ct1, dp0);
#
#         end % function Track
#
#
#         function TrackLibrary(drift,trackingMethod,libroutine)
#
#             drift1.length = drift.length;
#             if(~isempty(drift.aperture))
#                drift1.apertureX = drift.aperture(1);
#                drift1.apertureY = drift.aperture(2);
#             end
#
#             calllib(trackingMethod,[libroutine 'Drift'],drift1);
#
#         end % function TrackLibrary
#
#     end % methods
#
# end % classdef Drift
