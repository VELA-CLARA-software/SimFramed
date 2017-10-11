# classdef Sextupole < handle
#     % Sextupole class
#     %
#     % Properties:
#     %   name
#     %   length
#     %   gradient
#     %   aperture
#     %
#     % Methods:
#     %   Track
#     %   TrackSpin
#     %   GetBField
#
#     properties
#         name     = ''; % string
#         length   = 0;  % in metres
#         gradient = 0;  % in tesla/metre^2
#         aperture = []; % 1x2 array of elliptical aperture half-axes, in metres
#     end % properties
#
#     methods
#
#
#         function beam = Track(sextupole,beam)
#             % beam2 = Sextupole.Track(beam1)
#             % Applies the transfer map for a sextupole to the particles
#             % in beam1.
#
#             [x0, px0, y0, py0, ct0, dp0] = beam.GetParticles();
#
#             beta0  = beam.beta;
#
#             ds = sextupole.length;
#             k2 = sextupole.gradient / beam.rigidity; % normalised gradient
#
#             % First apply a drift through ds/2
#             d1  = sqrt(1 - px0.*px0 - py0.*py0 + 2*dp0/beta0 + dp0.*dp0);
#
#             x1  = x0  + ds*px0./d1/2;
#             y1  = y0  + ds*py0./d1/2;
#             ct1 = ct0 + ds*(1 - (1 + beta0*dp0)./d1)/beta0/2;
#
#             % Next, apply a sextupole 'kick'
#             px1 = px0 - (x1.*x1 - y1.*y1)*k2*ds/2;
#             py1 = py0 + x1.*y1*k2*ds;
#
#             % Finally, apply a second drift through ds/2
#             d1  = sqrt(1 - px1.*px1 - py1.*py1 + 2*dp0/beta0 + dp0.*dp0);
#
#             x2  = x1  + ds*px1./d1/2;
#             y2  = y1  + ds*py1./d1/2;
#             ct2 = ct1 + ds*(1 - (1 + beta0*dp0)./d1)/beta0/2;
#
#             beam.SetParticles(x2, px1, y2, py1, ct2, dp0);
#
#         end % function Track
#
#
#         function TrackLibrary(sextupole,trackingMethod,libroutine)
#
#             sextupole1.length   = sextupole.length;
#             sextupole1.gradient = sextupole.gradient;
#             if(~isempty(sextupole.aperture))
#                sextupole1.apertureX = sextupole.aperture(1);
#                sextupole1.apertureY = sextupole.aperture(2);
#             end
#
#             calllib(trackingMethod,[libroutine 'Sextupole'],sextupole1);
#
#         end % function TrackLibrary
#
#
#         function [bx, by, bz] = GetBField(sextupole,beam)
#             % [bx, by, bz] = Sextupole.GetBField(beam)
#             % Returns the magnetic field (in tesla) at the locations of
#             % the particles in the beam.
#
#             [x, ~, y] = beam.GetParticles();
#
#             bx = sextupole.gradient * (x.*x - y.*y) / 2;
#             by = sextupole.gradient * x .* y;
#             bz = 0;
#
#         end % function GetBField
#
#
#         function beam = TrackSpin(sextupole,beam)
#             % beam2 = Sextupole.Trackspin(beam1)
#             % Tracks particle spins through a Sextupole.
#
#             [bx, by, bz] = sextupole.GetBField(beam);
#
#             Utilities.SpinRotation(beam,bx,by,bz,sextupole.length);
#
#         end % function TrackSpin
#
#
#     end % methods
#
# end % classdef Sextupole
