# SAM to Python Conversion
# DJS August 2017
# Version 0.1
#
from ComponentBase import ComponentBase

class Dipole(ComponentBase):
    def __init__(self,field=0,gradient=0,hgap=0,e1=0,fint1=0,e2=0):
        # dipole field strength, in tesla
        self.field = field
        # quadrupole gradient, in tesla/metre
        self.gradient = gradient
        # gap height (for fringe field map), in metres
        self.hgap = hgap
        # entrance pole face rotation, in radians
        self.e1 = e1
        # entrance fringe field integral
        self.fint1 = fint1
        # exit pole face rotation, in radians
        self.e2 = e2
        # exit fringe field integral
        self.fint2 = fint2
        #  bending angle in radians
        self.__theta = 0

    # set the bend angle based on radius of curvature
    def set_curvature(self, curvature):
        self.__theta = self.length * curvature

    # get the radius of curvature
    def get_curvature(self):
        return self.__theta / self.length


# classdef Dipole < handle
#     % Dipole class
#     %
#     % Properties:
#     %   name
#     %   length
#     %   field
#     %   gradient
#     %   hgap
#     %   e1
#     %   fint1
#     %   e2
#     %   fint2
#     %   aperture
#     %
#     % Methods:
#     %   Track
#     %   TrackSpin
#     %   GetBField
#
#     properties
#         name      = ''; % string
#         length    = 0;  % dipole length, in metres
#         field     = 0;  % dipole field strength, in tesla
#         gradient  = 0;  % quadrupole gradient, in tesla/metre
#         hgap      = 0;  % gap height (for fringe field map), in metres
#         e1        = 0;  % entrance pole face rotation, in radians
#         fint1     = 0;  % entrance fringe field integral
#         e2        = 0;  % exit pole face rotation, in radians
#         fint2     = 0;  % exit fringe field integral
#         aperture  = []; % 1x2 array of elliptical aperture half-axes, in metres
#     end % properties
#
#     properties (Access=private)
#         theta     = 0; % bending angle in radians
#     end % properties (private)
#
#     properties (Dependent=true)
#         angle;         % bending angle in radians
#         curvature;     % reciprocal of radius of curvature, in /metre
#     end % properties (dependent)
#
#
#     methods
#
#         function set.angle(dipole,angle)
#             dipole.theta = angle;
#         end % set angle method
#
#         function angle = get.angle(dipole)
#             angle = dipole.theta;
#         end % get angle method
#
#
#         function set.curvature(dipole,curvature)
#             dipole.theta = dipole.length * curvature;
#         end % set curvature method
#
#         function curvature = get.curvature(dipole)
#             curvature = dipole.theta / dipole.length;
#         end % get curvature method
#
#
#         function beam = Track(dipole,beam)
#             % beam2 = Dipole.Track(beam1)
#             % Applies the transfer map for a dipole to the particles in
#             % in beam1.
#
#             [x0, px0, y0, py0, ct0, dp0] = beam.GetParticles();
#
#             beta0  = beam.beta;
#
#             ds     = dipole.length;
#             k0     = dipole.field / beam.rigidity;
#             d1     = sqrt(1 + 2*dp0/beta0 + dp0.*dp0);
#
#             % First, apply a map for the entrance fringe field
#             sine1 = sin(dipole.e1);
#             phi   = 2*dipole.fint1*dipole.hgap*k0*(1+sine1*sine1)/cos(dipole.e1);
#             r10   = k0*tan(dipole.e1);
#             r32   =-k0*tan(dipole.e1 - phi);
#
#             px1   = px0 + r10*x0;
#             py1   = py0 + r32*y0;
#
#             % Then, apply a map for the body of the dipole
#             h   = dipole.curvature;
#             k1  = dipole.gradient / beam.rigidity;
#             a1  = h - k0./d1;
#
#             wx  = sqrt((h*k0 + k1)./d1);
#             xc  = cos(wx*ds);
#             xs  = sin(wx*ds)./wx;
#             xs2 = sin(2*wx*ds)./wx;
#
#             wy  = sqrt(k1./d1);
#             yc  = cosh(wy*ds);
#             ys  = ds;
#             ys2 = 2*ds;
#
#             if(wy~=0)
#                 ys  = sinh(wy*ds)./wy;
#                 ys2 = sinh(2*wy*ds)./wy;
#             end
#
#             x2  =             x0.*xc + px1.*xs./d1 + a1.*(1-xc)./wx./wx;
#             px2 =-d1.*wx.*wx.*x0.*xs + px1.*xc     + a1.*xs.*d1;
#
#             y2  =             y0.*yc + py1.*ys./d1;
#             py2 = d1.*wy.*wy.*y0.*ys + py1.*yc;
#
#             d0  = 1/beta0 + dp0;
#
#             c0  = (1/beta0 - d0./d1)*ds - ...
#                   d0.*a1.*(h*(ds-xs) + ...
#                            a1.*(2*ds-xs2)/8)./wx./wx./d1;
#
#             c1  =-d0.*(h*xs - ...
#                        a1.*(2*ds-xs2)/4)./d1;
#
#             c2  =-d0.*(h*(1-xc)./wx./wx + ...
#                        a1.*xs.*xs/2)./d1./d1;
#
#             c11 =-d0.*wx.*wx.*(2*ds-xs2)./d1/8;
#             c12 = d0.*wx.*wx.*xs.*xs./d1./d1/2;
#             c22 =-d0.*(2*ds+xs2)./d1./d1./d1/8;
#
#             c33 =-d0.*wy.*wy.*(2*ds-ys2)./d1/8;
#             c34 =-d0.*wy.*wy.*ys.*ys./d1./d1/2;
#             c44 =-d0.*(2*ds+ys2)./d1./d1./d1/8;
#
#             ct1 = ct0 + c0 + ...
#                         c1.*x0      + c2.*px1 + ...
#                         c11.*x0.*x0 + c12.*x0.*px1 + c22.*px1.*px1 + ...
#                         c33.*y0.*y0 + c34.*y0.*py1 + c44.*py1.*py1;
#
#
#             % Finally, apply a map for the exit fringe field
#             sine2 = sin(dipole.e2);
#             phi = 2*dipole.fint2*dipole.hgap*k0*(1+sine2*sine2)/cos(dipole.e2);
#             r10 = k0*tan(dipole.e2);
#             r32 =-k0*tan(dipole.e2 - phi);
#
#             px3 = px2 + r10*x2;
#             py3 = py2 + r32*y2;
#
#             beam.SetParticles(x2, px3, y2, py3, ct1, dp0);
#
#         end % function Track
#
#
#         function TrackLibrary(dipole,trackingMethod,libroutine)
#
#             dipole1.length    = dipole.length;
#             dipole1.field     = dipole.field;
#             dipole1.gradient  = dipole.gradient;
#             dipole1.curvature = dipole.curvature;
#             dipole1.hgap      = dipole.hgap;
#             dipole1.e1        = dipole.e1;
#             dipole1.fint1     = dipole.fint1;
#             dipole1.e2        = dipole.e2;
#             dipole1.fint2     = dipole.fint2;
#
#             if(~isempty(dipole.aperture))
#                dipole1.apertureX = dipole.aperture(1);
#                dipole1.apertureY = dipole.aperture(2);
#             end
#
#             calllib(trackingMethod,[libroutine 'Dipole'],dipole1);
#
#         end % function TrackC
#
#
#         function [bx, by, bz] = GetBField(dipole,beam)
#             % [bx, by, bz] = Dipole.GetBField(beam)
#             % Returns the magnetic field (in tesla) at the locations of
#             % the particles in the beam.
#
#             [x, ~, y] = beam.GetParticles();
#
#             bx = y*dipole.gradient;
#             by = x*dipole.gradient + dipole.field;
#             bz = 0;
#
#         end % function GetBField
#
#
#         function beam = TrackSpin(dipole,beam)
#             % beam2 = Dipole.Trackspin(beam1)
#             % Tracks particle spins through a Dipole.
#
#             [bx, by, bz] = dipole.GetBField(beam);
#
#             Utilities.SpinRotation(beam,bx,by,bz,dipole.length);
#
#             % Account for the rotation of the local coordinate system
#
#             [theta1, phi1] = beam.GetSpins();
#
#             polsnx2    = sin(theta1).*cos(phi1);
#             polsny2    = sin(theta1).*sin(phi1);
#             polsnz2    = cos(theta1);
#
#             thetaprime = acos(polsny2);
#             phiprime   = atan2(polsnx2,polsnz2) + (dipole.length*dipole.curvature);
#
#             polsnx3    = sin(thetaprime).*sin(phiprime);
#             polsny3    = cos(thetaprime);
#             polsnz3    = sin(thetaprime).*cos(phiprime);
#
#             phi2       = atan2(polsny3,polsnx3);
#             theta2     = acos(polsnz3);
#
#             beam.SetSpins(theta2, phi2);
#
#         end % function TrackSpin
#
#
#     end % methods
#
# end % classdef Dipole
