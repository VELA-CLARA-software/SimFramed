# SAM to Python Conversion
# DJS August 2017
# Version 0.1
#
from ComponentBase import ComponentBase
from ..SAMMlab import Beam
import numpy

class MasterOscillator(ComponentBase):
    def __init__(self,voltage=0,harmonic=1,phase=numpy.pi,length=0, name="", aperture=[]):
        ComponentBase.__init__(self, length, name, aperture)
        # in volts
        self.volts = volts

# classdef MasterOscillator
#     % MasterOscillator (abstract) class
#     %
#     % Methods:
#     %   SetFrequency
#     %   GetFrequency
#
#     methods (Static)
#
#         function frequency = SetFrequency(frequency)
#             % Set the frequency of the master oscillator
#             global MasterOscillatorFrequency;
#             MasterOscillatorFrequency = frequency;
#         end
#
#         function frequency = GetFrequency()
#             % Get the frequency of the master oscillator
#             global MasterOscillatorFrequency;
#             frequency = MasterOscillatorFrequency;
#         end
#
#     end
#
# end