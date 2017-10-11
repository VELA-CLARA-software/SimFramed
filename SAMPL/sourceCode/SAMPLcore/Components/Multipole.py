# classdef Multipole < handle
#     % Multipole class
#     %
#     % Properties:
#     %   name
#     %   angle
#     %   curvature
#     %   field
#     %   aperture
#     %
#     % Methods:
#     %   Track
#     %   TrackSpin
#     %   GetBField
#
#     properties
#         name      = ''; % string
#         angle     = 0;  % bending angle of reference trajectory, in radians
#         curvature = 0;  % reciprocal of radius of curvature of reference trajectory, in 1/metres (for applying weak focusing)
#         field     = []; % Nx2 array of N field components [index normal-component+i*skew-component]
#         aperture  = []; % 1x2 array of elliptical aperture half-axes, in metres
#     end % properties
#
#     properties (Constant=true)
#         length   = 0;
#     end % properties (constant)
#
#     methods
#
#
#         function beam = Track(multipole,beam)
#             % beam2 = Multipole.Track(beam1)
#             % Applies the transfer map for a multipole to the particles
#             % in beam1.
#
#             [x0, px0, y0, py0, ct0, dp0] = beam.GetParticles();
#
#             dpx = -multipole.angle * (1 + dp0 - multipole.curvature * x0);
#
#             % Apply a multipole 'kick'
#             xi0 = x0 + 1i*y0;
#
#             for indx = 1:size(multipole.field,1)
#                 n   = multipole.field(indx,1);
#                 dpx = dpx + multipole.field(indx,2) * (xi0.^n) / beam.rigidity / factorial(n);
#             end
#
#             px1 = px0 - real(dpx);
#             py1 = py0 + imag(dpx);
#             ct1 = ct0 - multipole.angle * x0;
#
#             beam.SetParticles(x0, px1, y0, py1, ct1, dp0);
#
#         end % function Track
#
#
#         function TrackLibrary(multipole,trackingMethod,libroutine)
#
#             multipole.field       = sortrows(multipole.field,1);
#
#             multipole1.angle      = multipole.angle;
#             multipole1.curvature  = multipole.curvature;
#             multipole1.ncpts      = size(multipole.field,1);
#             multipole1.fieldindex = libpointer('doublePtr',real(multipole.field(:,1)));
#             multipole1.bnL        = libpointer('doublePtr',real(multipole.field(:,2)));
#             multipole1.anL        = libpointer('doublePtr',imag(multipole.field(:,2)));
#             if(~isempty(multipole.aperture))
#                multipole1.apertureX = multipole.aperture(1);
#                multipole1.apertureY = multipole.aperture(2);
#             end
#
#             calllib(trackingMethod,[libroutine 'Multipole'],multipole1);
#
#         end % function TrackLibrary
#
#
#         function [bx, by, bz] = GetBField(multipole,beam)
#             % [bx, by, bz] = Multipole.GetBField(beam)
#             % Returns the *integrated* magnetic field (in tesla metres) at
#             % the locations of the particles in the beam.
#
#             [x, ~, y] = beam.GetParticles();
#
#             xi = x + 1i * y;
#             b  = 0;
#
#             for indx = 1:size(multipole.field,1)
#                 n = multipole.field(indx,1);
#                 b = b + multipole.field(indx,2) * xi.^n / factorial(n);
#             end
#
#             bx = real(b);
#             by = imag(b);
#             bz = 0;
#
#         end % function GetBField
#
#
#         function beam = TrackSpin(multipole,beam)
#             % beam2 = Multipole.Trackspin(beam1)
#             % Tracks particle spins through a Multipole.
#
#             [bxL, byL, bzL] = multipole.GetBField(beam);
#
#             nomlength = 1e-9;
#
#             Utilities.SpinRotation(beam, bxL/nomlength, byL/nomlength, bzL/nomlength, nomlength);
#
#             % Account for the rotation of the local coordinate system
#
#             if multipole.angle ~= 0
#
#                 [theta1, phi1] = beam.GetSpins();
#
#                 polsnx2    = sin(theta1).*cos(phi1);
#                 polsny2    = sin(theta1).*sin(phi1);
#                 polsnz2    = cos(theta1);
#
#                 thetaprime = acos(polsny2);
#                 phiprime   = atan2(polsnx2,polsnz2) + multipole.angle;
#
#                 polsnx3    = sin(thetaprime).*sin(phiprime);
#                 polsny3    = cos(thetaprime);
#                 polsnz3    = sin(thetaprime).*cos(phiprime);
#
#                 phi2       = atan2(polsny3,polsnx3);
#                 theta2     = acos(polsnz3);
#
#                 beam.SetSpins(theta2, phi2);
#
#             end
#
#         end % function TrackSpin
#
#
#     end % methods
#
# end % classdef Multipole
