# classdef Marker < handle
#     % Marker class
#     %
#     % Properties:
#     %   name
#     %
#     % Methods:
#     %   Track
#
#     properties
#         name     = ''; % string
#     end % properties
#
#     properties (Constant)
#         length   = 0;  % in metres
#         aperture = []; % 1x2 array of elliptical aperture half-axes, in metres
#     end % constant properties
#
#     methods
#
#
#         function beam = Track(~,beam)
#             % beam2 = Marker.Track(beam1)
#             % Applies the transfer map for a marker to the particles
#             % in beam1.
#
#         end % function Track
#
#
#         function TrackLibrary(~,trackingMethod,libroutine)
#
#             calllib(trackingMethod,[libroutine 'Marker']);
#
#         end % function TrackLibrary
#
#
#     end % methods
#
#
# end % classdef Marker
