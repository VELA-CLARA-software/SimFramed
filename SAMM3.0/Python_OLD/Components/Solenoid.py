# classdef Solenoid < handle
#     % Solenoid class
#     %
#     % Properties:
#     %   name
#     %   length
#     %   field
#     %   taper
#     %   aperture
#     %
#     % Methods:
#     %   Track
#     %   TrackSpin
#     %   GetBField
#
#     properties
#         name     = ''; % string
#         length   = 0;  % solenoid length in metres
#         field    = 0;  % solenoid field in tesla
#         taper    = 0;  % taper partemeter, in /metre
#         aperture = []; % 1x2 array of elliptical aperture half-axes, in metres
#     end % properties
#
#     methods
#
#         function beam = Track(solenoid,beam)
#             % beam2 = Solenoid.Track(beam1)
#             % Applies the dynamical map for a solenoid to the particles
#             % in beam1.
#             % Note that the path length is approximate.
#
#             [x0, px0, y0, py0, ct0, dp0] = beam.GetParticles();
#
#             ds    = solenoid.length;
#             b0    = solenoid.field;
#             g     = solenoid.taper;
#
#             brho  = beam.rigidity;
#
#             beta0 = beam.beta;
#
#             % First, apply a drift map through ds/2
#             % to the longitudinal coordinate
#             d1    = sqrt(1 - px0.*px0 - py0.*py0 + 2*dp0/beta0 + dp0.*dp0);
#             ct1   = ct0 + ds*(1 - (1 + beta0*dp0)./d1)/beta0/2;
#
#             % Now apply a solenoid map to the transverse variables
#             gds   = g*ds;
#             b1    = b0 / (1 + gds);
#
#             helmA = sqrt(b0/b1) * ones(size(x0));
#             helmB = 2*brho*(1+dp0)/sqrt(b0*b1);
#             helmF = -1 ./ helmB;
#             helmG =  1 ./ helmA;
#
#             w = b0*ds/2/brho./(1+dp0);
#             if gds~=0
#                 w = w*log(1+gds)/gds;
#             end;
#
#             cw2 = cos(w).*cos(w);
#             s2w = sin(2*w);
#             sw2 = sin(w).*sin(w);
#
#             x1  = helmA.*cw2.*x0   + helmB.*s2w.*px0/2 + helmA.*s2w.*y0/2 + helmB.*sw2.*py0;
#             px1 = helmF.*s2w.*x0/2 + helmG.*cw2.*px0   + helmF.*sw2.*y0   + helmG.*s2w.*py0/2;
#             y1  =-helmA.*s2w.*x0/2 - helmB.*sw2.*px0   + helmA.*cw2.*y0   + helmB.*s2w.*py0/2;
#             py1 =-helmF.*sw2.*x0   - helmG.*s2w.*px0/2 + helmF.*s2w.*y0/2 + helmG.*cw2.*py0;
#
#             % First, apply a drift map through ds/2
#             % to the longitudinal coordinate
#             d1  = sqrt(1 - px1.*px1 - py1.*py1 + 2*dp0/beta0 + dp0.*dp0);
#             ct2 = ct1 + ds*(1 - (1 + beta0*dp0)./d1)/beta0/2;
#
#             % Set the new values for the dynamical variables
#             beam.SetParticles(x1, px1, y1, py1, ct2, dp0);
#
#         end % function Track
#
#
#         function TrackLibrary(solenoid,trackingMethod,libroutine)
#
#             solenoid1.length = solenoid.length;
#             solenoid1.taper  = solenoid.taper;
#             solenoid1.field  = solenoid.field;
#             if(~isempty(solenoid.aperture))
#                solenoid1.apertureX = solenoid.aperture(1);
#                solenoid1.apertureY = solenoid.aperture(2);
#             end
#
#             calllib(trackingMethod,[libroutine 'Solenoid'],solenoid1);
#
#         end % function TrackLibrary
#
#
#         function [bx, by, bz] = GetBField(solenoid,beam)
#             % [bx, by, bz] = Solenoid.GetBField(beam)
#             % Returns the magnetic field (in tesla) at the locations of
#             % the particles in the beam.
#
#             [x, ~, y] = beam.GetParticles();
#
#             ds  = solenoid.length;
#             b0  = solenoid.field;
#             g   = solenoid.taper;
#
#             gds = g*ds;
#             bx  = x*b0*g/(1+gds)/2;
#             by  = y*b0*g/(1+gds)/2;
#             bz  = b0*ones(size(x));
#
#             if gds~=0
#                 bz = bz*log(1+gds)/gds;
#             end
#
#         end % function GetBField
#
#
#         function beam = TrackSpin(solenoid,beam)
#             % beam2 = Solenoid.Trackspin(beam1)
#             % Tracks particle spins through a Solenoid.
#
#             [bx, by, bz] = solenoid.GetBField(beam);
#
#             Utilities.SpinRotation(beam,bx,by,bz,solenoid.length);
#
#         end % function TrackSpin
#
#     end % methods
#
# end % classdef Solenoid
