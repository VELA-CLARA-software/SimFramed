# SAM to Python Conversion
# DJS August 2017
# Version 0.1
#
import ../SAMMlab/PhysicalConstants

class Electron(object):
    def __init__(self,componentlist=[], tracklocal=False):
        self.__charge = -PhysicalConstants.PositronCharge
        self.__mass   =  PhysicalConstants.ElectronMass
        self.__Cq     =  55*PhysicalConstants.PlanckConstant/32/sqrt(3)/2/pi/PhysicalConstants.SpeedOfLight/PhysicalConstants.ElectronMass

        self.__Cgamma = PhysicalConstants.PositronCharge^2/3/\
                        PhysicalConstants.VacuumPermittivity/\
                        (PhysicalConstants.ElectronMass*PhysicalConstants.SpeedOfLight^2)^4;

# classdef Electron
#     % Electron
#
#     properties (Constant)
#
#         charge = -PhysicalConstants.PositronCharge;
#         mass   =  PhysicalConstants.ElectronMass;
#         g      =  PhysicalConstants.PositrongFactor;
#
#         Cq     = 55*PhysicalConstants.PlanckConstant/32/sqrt(3)/2/pi/PhysicalConstants.SpeedOfLight/PhysicalConstants.ElectronMass;
#         Cgamma = PhysicalConstants.PositronCharge^2/3/PhysicalConstants.VacuumPermittivity/(PhysicalConstants.ElectronMass*PhysicalConstants.SpeedOfLight^2)^4;
#
#     end % properties
#
# end % classdef