

# classdef FieldMap < handle
#     % FieldMap class
#     %
#     % Properties:
#     %   name
#     %   length
#     %   aperture
#     %
#     % Methods:
#     %   Track
#
#     properties
#         name     = ''; % string
#         length   = 0;  % in metres
#         aperture = []; % 1x2 array of elliptical aperture half-axes, in metres
#         gridX    = []; % x coordinates of grid points
#         gridY    = []; % y coordinates of grid points
#         gridZ    = []; % z coordinates of grid points
#         Bx       = []; % horizontal field component at the grid points
#         By       = []; % vertical field component at the grid points
#         Bz       = []; % longitudinal field component at the grid points
#         interpmethod = 'linear'; % interpolation method for field (see Matlab help interp3 for options)
#         nsteps   = 10; % number of steps to take in tracking through the field
#     end % properties
#
#     methods
#
#         function beam = Track(fieldmap,beam)
#             % beam2 = FieldMap.Track(beam1)
#             % Applies the transfer map for a drift to the particles
#             % in beam1.
#
#             [x0, px0, y0, py0, ct0, dp0] = beam.GetParticles();
#
#             z0      = zeros(size(x0));
#
#             species = beam.species;
#
#             P0      = beam.momentum/species.mass/PhysicalConstants.SpeedOfLight;
#             beta0   = beam.beta;
#             gamma0  = beam.gamma;
#
#             gamma   = gamma0*(1 + beta0*dp0);
#
#             b0      = sqrt(1 - gamma.^(-2));
#             bx0     = P0*px0./gamma;
#             by0     = P0*py0./gamma;
#             bz0     = sqrt(b0.*b0 - bx0.*bx0 - by0.*by0);
#
#             k       = (species.charge/species.mass/PhysicalConstants.SpeedOfLight)./gamma;
#
#             for n = 1:fieldmap.nsteps
#
#                 cdt   = (fieldmap.length - z0) .* ones(size(x0)) / (fieldmap.nsteps + 1 - n) ./ bz0;
#
#                 % Note the order of the coordinates (y,x,z) when calling interp3
#                 % This is a Matlab peculiarity!
#                 Bx0   = k.*interp3(fieldmap.gridY,fieldmap.gridX,fieldmap.gridZ,fieldmap.Bx,y0,x0,z0,fieldmap.interpmethod);
#                 By0   = k.*interp3(fieldmap.gridY,fieldmap.gridX,fieldmap.gridZ,fieldmap.By,y0,x0,z0,fieldmap.interpmethod);
#                 Bz0   = k.*interp3(fieldmap.gridY,fieldmap.gridX,fieldmap.gridZ,fieldmap.Bz,y0,x0,z0,fieldmap.interpmethod);
#
#                 bmag  = sqrt(Bx0.*Bx0 + By0.*By0 + Bz0.*Bz0);
#
#                 bdotv = Bx0.*bx0 + By0.*by0 + Bz0.*bz0;
#                 s1    = sin(bmag.*cdt)./bmag;
#                 s2    = bdotv.*(bmag.*cdt - sin(bmag.*cdt)).*bmag.^(-3);
#                 c1    = cos(bmag.*cdt);
#                 c2    = (1 - c1).*bmag.^(-2);
#
#                 x1  = x0 + bx0.*s1 + (by0.*Bz0 - bz0.*By0).*c2 + Bx0.*s2;
#                 y1  = y0 + by0.*s1 + (bz0.*Bx0 - bx0.*Bz0).*c2 + By0.*s2;
#                 z1  = z0 + bz0.*s1 + (bx0.*By0 - by0.*Bx0).*c2 + Bz0.*s2;
#                 ct0 = ct0 + (bz0/beta0 - 1) .* cdt;
#
#                 bx1 = bx0.*c1 + (by0.*Bz0 - bz0.*By0).*s1 + Bx0.*bdotv.*c2;
#                 by1 = by0.*c1 + (bz0.*Bx0 - bx0.*Bz0).*s1 + By0.*bdotv.*c2;
#                 bz1 = bz0.*c1 + (bx0.*By0 - by0.*Bx0).*s1 + Bz0.*bdotv.*c2;
#
#                 x0  = x1;
#                 y0  = y1;
#                 z0  = z1;
#
#                 bx0 = bx1;
#                 by0 = by1;
#                 bz0 = bz1;
#
#             end
#
#             px0 = bx0.*gamma/P0;
#             py0 = by0.*gamma/P0;
#
#             beam.SetParticles(x0, px0, y0, py0, ct0, dp0);
#
#         end % function Track
#
#
#         function TrackLibrary(fieldmap,trackingMethod,libroutine)
#
#             fieldmap1.length = fieldmap.length;
#             if(~isempty(fieldmap.aperture))
#                fieldmap1.apertureX = fieldmap.aperture(1);
#                fieldmap1.apertureY = fieldmap.aperture(2);
#             end
#
#             fieldmap1.gridXsize = size(fieldmap.gridX,2);
#             fieldmap1.gridYsize = size(fieldmap.gridY,2);
#             fieldmap1.gridZsize = size(fieldmap.gridZ,2);
#
#             fieldmap1.gridX = libpointer('doublePtr',fieldmap.gridX);
#             fieldmap1.gridY = libpointer('doublePtr',fieldmap.gridY);
#             fieldmap1.gridZ = libpointer('doublePtr',fieldmap.gridZ);
#
#             fieldmap1.Bx = libpointer('doublePtr',fieldmap.Bx(1:end));
#             fieldmap1.By = libpointer('doublePtr',fieldmap.By(1:end));
#             fieldmap1.Bz = libpointer('doublePtr',fieldmap.Bz(1:end));
#
#             fieldmap1.nsteps = fieldmap.nsteps;
#
#             calllib(trackingMethod,[libroutine 'FieldMap'],fieldmap1);
#
#         end % function TrackLibrary
#
#     end % methods
#
# end % classdef FieldMap
