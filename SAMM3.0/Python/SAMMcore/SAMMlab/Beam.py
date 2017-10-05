# SAM to Python Conversion
# DJS August 2017
# Version 0.1
#
import math
from ..Particles import Electron

from ..SAMMlab import PhysicalConstants
from ..SAMMlab import PhysicalUnits
import numpy as np
import math


class Beam(object):
    def __init__(self, species = Electron.Electron, bunchcharge = 0, distance = [],
                 globaltime = 0,  spins = 0, energy = 0,
                 momentum = 0,  rigidity = 0):
        # particle type, e.g. electron, positron, proton
        self.species = species()
        # bunch charge in coulombs
        self.bunchcharge = bunchcharge
        # 2xN array, recording total distance travelled by each particle%
        # this is s - the independant coordinate
        self.distance = np.array(distance)
        # global time in seconds
        self.globaltime = globaltime
        # 2xN array, polar angles of spin vector of N particles
        self.spins = spins
        # 6xN array, phase space variables of N particles
        #self.particles = particles

        # reference energy in joules
        if energy != 0:
            print 'calling energy setter??'
            self.energy = energy

        # reference momentum in kilogram metres/second
        if momentum != 0:
            self.momentum = momentum

        # beam rigidity in tesla metres
        # this is called last - if energy and momentum are also passed and not
        # consistent with rigidity, rigidity is assumed to be correct
        if rigidity != 0:
            self.__rigidity = rigidity

        # do we really need this... ?
        self.__brho = 0  # ??
        # Dependent=true properties get a __ prefix
        # relativistic beta (v/c)
        self.__beta = 0
        # relativistic gamma factor
        self.__gamma = 0
        # beam rigidity in tesla metres
        #self.__brho = 0
        # phase space coordinates of particles
        self.__phasespace = 0
        # why not have phase space as explicit lists? that are accessed anywhere?
        self.x  = np.array([1,1])
        self.px = np.array([2,2])
        self.y  = np.array([3,3])
        self.py = np.array([4,4])
        self.ct = np.array([5,5])
        self.dp = np.array([6,7])
        #self.__particles = np.array([self.x, self.px, self.y, self.py, self.ct, self.dp])
        self.__particles = np.column_stack((self.x, self.px, self.y, self.py, self.ct, self.dp))
        # for i in range(len(self.__particles)):
        #     print self.__particles[i]

#https://softwareengineering.stackexchange.com/questions/283514/python-best-way-to-have-interdependant-variable-in-a-class

    @property
    def particles(self):
        return np.transpose( [self.x, self.px, self.y, self.py, self.ct, self.dp] )

    @particles.setter
    def particles(self, particles ):
        if type(particles).__module__ == np.__name__:
            #print 'setting particles'
            [self.x, self.px, self.y, self.py, self.ct, self.dp] = np.column_stack((particles))
            #print self.x
            #print self.px
            #print self.y
            #print self.py
            # distance is whether they are in the aperture or not at each stage (!)
            self.distance = np.array([[0.0] * len(particles), [1.0] * len(particles)])
            #print self.distance
        else:
            print 'particles should be a numpy array'

    def gamma(self):
        return self.__gamma

    @property
    def energy(self):
        return math.sqrt(math.pow(self.__rigidity * self.species.charge, 2) * PhysicalConstants.SpeedOfLight2
                         + (self.species.mass2 * PhysicalConstants.SpeedOfLight4))

    @energy.setter
    def energy(self, energy):
        self.__rigidity = math.sqrt(energy**2 / PhysicalConstants.SpeedOfLight2
        - (self.species.mass2 * PhysicalConstants.SpeedOfLight2)) / self.species.charge
        # print('energy: set rigidity = ', self.__rigidity)

    @property
    def momentum(self):
        return  self.__rigidity * self.species.charge

    @momentum.setter
    def momentum(self, momentum):
        self.__rigidity = momentum / self.species.charge

    @property
    def rigidity(self):
        return self.__rigidity

    @rigidity.setter
    def rigidity(self,rigidity):
        self.__rigidity = abs(rigidity) * numpy.sign(self.species.charge)

        #if ~isempty(beam.phasespace) && beam.brho~=0
        if not self.phasespace and self.__brho != 0:

            # if we have some particles and we are scaling their momentum
            # (maybe after emitting radiation, or passing though an
            # rf structure? )


            P0 = self.__brho * self.species.charge;
            P1 = self.__rigidity  * self.species.charge

            E0 = math.sqrt(math.pow(P0,2)* PhysicalConstants.SpeedOfLight2 \
             + self.species.mass2 * PhysicalConstants.SpeedOfLight4 )

            E1 = math.sqrt(math.pow(P1,2)* PhysicalConstants.SpeedOfLight2 \
             + self.species.mass2 * PhysicalConstants.SpeedOfLight4 )

            b0 = P0*PhysicalConstants.SpeedOfLight/E0
            b1 = P1*PhysicalConstants.SpeedOfLight/E1

            #self.phasespace([2,4],:) = self.phasespace([2,4],:)*P0/P1;
            #self.phasespace(6,:) = (self.phasespace(6,:) + 1/b0)*P0/P1 - 1/b1;

            #self.brho = self.__rigidity

    @property
    def beta(self):
        bg = self.__rigidity * self.species.charge / self.species.mass / PhysicalConstants.SpeedOfLight
        self.__beta = bg/ math.sqrt(1+bg*bg)
        return self.__beta

    @property
    def gamma(self):
        bg = self.__rigidity * self.species.charge / self.species.mass / PhysicalConstants.SpeedOfLight
        self.__gamma = math.sqrt(1+bg*bg)
        return self.__gamma



#
# classdef Beam < handle
#     % Beam class
#     %
#     % Properties:
#     %   species
#     %   bunchcharge
#     %   distance
#     %   globaltime
#     %   particles
#     %   spins
#     %   energy
#     %   momentum
#     %   rigidity
#     %   beta ('get' only)
#     %   gamma ('get' only)
#     %
#     % Methods:
#     %   SetParticles
#     %   GetParticles
#     %   SetSpins
#     %   GetSpins
#
#     properties
#         species;             % particle type, e.g. electron, positron, proton
#         bunchcharge = 0;     % bunch charge in coulombs
#         distance    = [];    % 2xN array, recording total distance travelled by each particle
#         globaltime  = 0;     % global time in seconds
#         spins       = [];    % 2xN array, polar angles of spin vector of N particles
#     end % properties
#
#     properties (Access=private)
#         phasespace;          % phase space coordinates of particles
#         brho;                % beam rigidity in tesla metres
#     end % properties (private)
#
#     properties (Dependent=true)
#         particles;           % 6xN array, phase space variables of N particles
#         energy;              % reference energy in joules
#         momentum;            % reference momentum in kilogram metres/second
#         rigidity;            % beam rigidity in tesla metres
#         beta;                % relativistic beta (v/c)
#         gamma;               % relativistic gamma factor
#     end % properties (dependent)
#
#
#     methods
#
#         function beam = Beam(species)
#             % initialise a beam with specified particle species
#             beam.species = species;
#         end % function Beam
#
#
#         function set.energy(beam,energy)
#             spec = beam.species;
#             beam.rigidity = sqrt(energy^2/PhysicalConstants.SpeedOfLight^2 ...
#                             - spec.mass^2*PhysicalConstants.SpeedOfLight^2)/spec.charge;
#         end % energy set method
#
#         function energy = get.energy(beam)
#             spec = beam.species;
#             energy = sqrt((beam.brho * spec.charge)^2*PhysicalConstants.SpeedOfLight^2 ...
#                                         + spec.mass^2*PhysicalConstants.SpeedOfLight^4);
#         end % energy get method
#
#
#         function set.momentum(beam,momentum)
#             spec = beam.species;
#             beam.rigidity = momentum/spec.charge;
#         end % momentum set method
#
#         function momentum = get.momentum(beam)
#             spec = beam.species;
#             momentum = beam.brho * spec.charge;
#         end % momentum get method
#
#
#         function set.rigidity(beam,rigidity)
#
#             spec = beam.species;
#             rigidity = abs(rigidity)*sign(spec.charge);
#
#             if ~isempty(beam.phasespace) && beam.brho~=0
#
#                 P0 = beam.brho * spec.charge;
#                 P1 = rigidity  * spec.charge;
#
#                 E0 = sqrt(P0^2*PhysicalConstants.SpeedOfLight^2 ...
#                  + spec.mass^2*PhysicalConstants.SpeedOfLight^4);
#                 E1 = sqrt(P1^2*PhysicalConstants.SpeedOfLight^2 ...
#                  + spec.mass^2*PhysicalConstants.SpeedOfLight^4);
#
#                 b0 = P0*PhysicalConstants.SpeedOfLight/E0;
#                 b1 = P1*PhysicalConstants.SpeedOfLight/E1;
#
#                 beam.phasespace([2,4],:) = beam.phasespace([2,4],:)*P0/P1;
#
#                 beam.phasespace(6,:) = (beam.phasespace(6,:) + 1/b0)*P0/P1 - 1/b1;
#
#             end
#
#             beam.brho = rigidity;
#
#         end % rigidity set method
#
#         function rigidity = get.rigidity(beam)
#             rigidity = beam.brho;
#         end % rigidity get method
#
#
#         function beta = get.beta(beam)
#             spec = beam.species;
#             bg = beam.brho * spec.charge / spec.mass / PhysicalConstants.SpeedOfLight;
#             beta = bg/sqrt(1+bg*bg);
#         end % beta get method
#
#
#         function gamma = get.gamma(beam)
#             spec = beam.species;
#             bg = beam.brho * spec.charge / spec.mass / PhysicalConstants.SpeedOfLight;
#             gamma = sqrt(1+bg*bg);
#         end % gamma get method
#
#
#         function set.particles(beam,particles)
#             beam.phasespace = particles;
#             beam.distance   = [ zeros(1,size(particles,2));...
#                                 ones( 1,size(particles,2)) ];
#         end % particles set method
#
#         function ps = get.particles(beam)
#             ps = beam.phasespace;
#         end % particles get method
#
#
#         function val = SetParticles(beam,x,px,y,py,ct,dp)
#             % Beam.SetParticles(x,px,y,py,ct,dp)
#
#             beam.phasespace(:,beam.distance(2,:)~=0) = [x; px; y; py; ct; dp];
#             val = 0;
#         end % SetParticles method
#
#         function [x, px, y, py, ct, dp] = GetParticles(beam)
#             % [x, px, y, py, ct, dp] = Beam.GetParticles()
#
#             inbeam = find(beam.distance(2,:));
#             x      = beam.phasespace(1,inbeam);
#             px     = beam.phasespace(2,inbeam);
#             y      = beam.phasespace(3,inbeam);
#             py     = beam.phasespace(4,inbeam);
#             ct     = beam.phasespace(5,inbeam);
#             dp     = beam.phasespace(6,inbeam);
#         end % GetParticles method
#
#
#         function val = SetSpins(beam,theta,phi)
#             % Beam.SetSpins(theta,phi)
#
#             beam.spins(:,beam.distance(2,:)~=0) = [theta; phi];
#             val = 0;
#         end % SetSpins method
#
#         function [theta, phi] = GetSpins(beam)
#             % [theta, phi] = Beam.GetSpins()
#
#             inbeam = find(beam.distance(2,:)~=0);
#             theta  = beam.spins(1,inbeam);
#             phi    = beam.spins(2,inbeam);
#         end % GetSpins method
#
#     end % methods
#
#
# end % classdef Beam
