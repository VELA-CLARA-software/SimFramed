# SAM to Python Conversion
# DJS August 2017
# Version 0.1
#
from ComponentBase import ComponentBase
# from ..SAMMlab import Beam
from ..SAMMlab import MasterOscillator
from ..SAMMlab import PhysicalConstants as PhyC
import numpy
import math

class RFCavity(ComponentBase):
    def __init__(self, voltage=0, harmonic=1, phase=numpy.pi, length=0,
                 name="", aperture=[]):
        ComponentBase.__init__(self, length, name, aperture)
        # in volts
        self.voltage = voltage
        # frequency (in units of master oscillator frequency)
        self.harmonic = harmonic
        # relative to master oscillator, in radians
        self.phase = phase
        # 1x2 array of elliptical aperture half-axes, in metres
        self.aperture = aperture

    # @property

    def setFrequency(self, f):
        MasterOscillator.frequency = f / self.harmonic

    def getFrequency(self):
        return self.harmonic * MasterOscillator.frequency

    def Track(self, beam):
        # Applies the dynamical map for an RF cavity to the particles
        # in beam1.  The reference momentum is unchanged.
        beta0 = beam.beta
        mofreq = MasterOscillator.frequency
        ds = self.length
        f = self.harmonic * mofreq
        # First apply a drift through ds/2
        d1 = math.sqrt(1 - beam.px**2 - beam.py**2 + 2 * beam.dp / beta0 +
                       beam.dp**2)
        x1 = beam.x + ds * beam.px / d1 / 2
        y1 = beam.y + ds * beam.py / d1 / 2
        ct1 = beam.ct + ds * (1 - (1 + beam.dp * beta0) / d1) / beta0 / 2
        # Next, apply an rf 'kick'
        p = math.floor(beam.globaltime * mofreq)
        beam.globaltime = beam.globaltime - (p / mofreq)
        t = beam.globaltime - (ct1 / (beta0 * PhyC.SpeedOfLight))
        ft = (f * t) - math.floor(f * t)
        vnorm = self.voltage / beam.rigidity / PhyC.SpeedOfLight
        dp1 = beam.dp + (vnorm * math.sin((2 * math.pi * ft) + self.phase))
        # Finally, apply a second drift through ds/2
        d1 = math.sqrt(1 - beam.px**2 - beam.py**2 + 2 * dp1 / beta0 +
                       dp1**2)
        beam.x = x1 + (ds * beam.x) / d1 / 2
        beam.y = y1 + (ds * beam.y) / d1 / 2
        beam.ct = ct1 + ds * (1 - (1 + beta0 * dp1) / d1) / beta0 / 2

        # save
        self.lastTrackedBeam = beam
    # @momentum.setter
    # def momentum(self, momentum):
    #    self.__rigidity = momentum / self.species.charge

    # function set.frequency(rfcavity,f)
    #        f1 = rfcavity.harmonic * MasterOscillator.GetFrequency();
    #        MasterOscillator.SetFrequency(MasterOscillator.GetFrequency()*f/f1);
    #    end

    #    function f = get.frequency(rfcavity)
    #        f = rfcavity.harmonic * MasterOscillator.GetFrequency();
    #    end


# classdef RFCavity < handle
#     % RFCavity class
#     %
#     % Properties:
#     %   name
#     %   length
#     %   voltage
#     %   harmonic
#     %   phase
#     %   aperture
#     %
#     % Methods:
#     %   Track
#
#     properties
#         name     = ''; % string
#         length   = 0;  % in metres
#         voltage  = 0;  % in volts
#         harmonic = 1;  % frequency (in units of master oscillator frequency)
#         phase    = pi; % relative to master oscillator, in radians
#         aperture = []; % 1x2 array of elliptical aperture half-axes, in metres
#     end % properties
#
#     properties (Dependent=true)
#         frequency;
#     end % properties (dependent)
#
#     methods
#
#         function set.frequency(rfcavity,f)
#             f1 = rfcavity.harmonic * MasterOscillator.GetFrequency();
#             MasterOscillator.SetFrequency(MasterOscillator.GetFrequency()*f/f1);
#         end
#
#         function f = get.frequency(rfcavity)
#             f = rfcavity.harmonic * MasterOscillator.GetFrequency();
#         end
#
#         function beam = Track(rfcavity,beam)
#             % beam2 = RFCavity.Track(beam1)
#             % Applies the dynamical map for an RF cavity to the particles
#             % in beam1.  The reference momentum is unchanged.
#
#             [x0, px0, y0, py0, ct0, dp0] = beam.GetParticles();
#
#             beta0 = beam.beta;
#
#             mofreq = MasterOscillator.GetFrequency();
#
#             ds  = rfcavity.length;
#             f   = rfcavity.harmonic*mofreq;
#
#             % First apply a drift through ds/2
#             d1  = sqrt(1 - px0.*px0 - py0.*py0 + 2*dp0/beta0 + dp0.*dp0);
#
#             x1  = x0  + ds*px0./d1/2;
#             y1  = y0  + ds*py0./d1/2;
#             ct1 = ct0 + ds*(1 - (1 + beta0*dp0)./d1)/beta0/2;
#
#             % Next, apply an rf 'kick'
#             p = floor(beam.globaltime*mofreq);
#             beam.globaltime = beam.globaltime - p/mofreq;
#             t = beam.globaltime - ct1/(beta0*PhysicalConstants.SpeedOfLight);
#             ft = f*t - floor(f*t);
#
#             vnorm = rfcavity.voltage/beam.rigidity/PhysicalConstants.SpeedOfLight;
#             dp1   = dp0 + vnorm*sin(2*pi*ft + rfcavity.phase);
#
#             % Finally, apply a second drift through ds/2
#             d1  = sqrt(1 - px0.*px0 - py0.*py0 + 2*dp1/beta0 + dp1.*dp1);
#
#             x2  = x1  + ds*px0./d1/2;
#             y2  = y1  + ds*py0./d1/2;
#             ct2 = ct1 + ds*(1 - (1 + beta0*dp1)./d1)/beta0/2;
#
#             beam.SetParticles(x2, px0, y2, py0, ct2, dp1);
#
#         end % function Track
#
#         function TrackLibrary(rfcavity,trackingMethod,libroutine)
#
#             mofreq = MasterOscillator.GetFrequency();
#
#             rfcavity1.length    = rfcavity.length;
#             rfcavity1.voltage   = rfcavity.voltage;
#             rfcavity1.frequency = rfcavity.harmonic * mofreq;
#             rfcavity1.phase     = rfcavity.phase;
#             rfcavity1.masteroscillatorfrequency = mofreq;
#
#             if(~isempty(rfcavity.aperture))
#                rfcavity1.apertureX = rfcavity.aperture(1);
#                rfcavity1.apertureY = rfcavity.aperture(2);
#             end
#
#             calllib(trackingMethod,[libroutine 'RFCavity'],rfcavity1);
#
#         end % function TrackLibrary
#
#     end % methods
#
# end % classdef RFCavity
