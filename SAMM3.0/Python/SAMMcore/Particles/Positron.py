# SAM to Python Conversion
# DJS August 2017
# Version 0.1
#
#from Particle import Particle
import math
from ..SAMMlab import PhysicalConstants


class Positron(object):
    def __init__(self):
        #Particle.__init__(self)
        #print 'Positron created'
        self.__charge = PhysicalConstants.PositronCharge
        self.__mass   = PhysicalConstants.ElectronMass
        self.__mass2  = self.__mass * self.__mass
        self.__Cq     = 55 * PhysicalConstants.PlanckConstant / 32 / math.sqrt(
            3) / 2 / math.pi / PhysicalConstants.SpeedOfLight / PhysicalConstants.ElectronMass

        self.__Cgamma = PhysicalConstants.PositronCharge2 / 3 / \
                        PhysicalConstants.VacuumPermittivity / \
                        math.pow( PhysicalConstants.ElectronMass * PhysicalConstants.SpeedOfLight2 , 4)

    @property
    def Positron(self):
        return Positron()

    @property
    def charge(self):
        return self.__charge

    @property
    def mass(self):
        return self.__mass

    @property
    def mass2(self):
        return self.__mass2

# classdef Positron
#     % Positron
#
#     properties (Constant)
#
#         charge =  PhysicalConstants.PositronCharge;
#         mass   =  PhysicalConstants.ElectronMass;
#         g      =  PhysicalConstants.PositrongFactor;
#
#         Cq     = 55*PhysicalConstants.PlanckConstant/32/sqrt(3)/2/pi/PhysicalConstants.SpeedOfLight/PhysicalConstants.ElectronMass;
#         Cgamma = PhysicalConstants.PositronCharge^2/3/PhysicalConstants.VacuumPermittivity/(PhysicalConstants.ElectronMass*PhysicalConstants.SpeedOfLight^2)^4;
#
#     end % properties
#
# end % classdef