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