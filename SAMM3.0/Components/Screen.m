classdef Screen < handle
    % Screen class
    %
    % Properties:
    %   name
    %   length
    %   aperture
    %
    % Methods:
    %   Track
% TP: I added this in because I want a screen component to SAMPLE so a had a fake one in SAMM
    properties
        name     = ''; % string
        length   = 0;  % in metres
        aperture = []; % 1x2 array of elliptical aperture half-axes, in metres
    end % properties

    methods

        function beam = Track(screen,beam)
            % beam2 = Screen.Track(beam1)
            % Applies the transfer map for a drift to the particles
            % in beam1.

            [x0, px0, y0, py0, ct0, dp0] = beam.GetParticles();
            % does nothing, not even calc x,y,sigmaX or sigmaY
            beam.SetParticles(x0, px0, y0, py0, ct0, dp0);

        end % function Track

        % the below wont exist
        function TrackLibrary(screen,trackingMethod,libroutine)

            screen1.length = screen.length;
            if(~isempty(screen.aperture))
               drift1.apertureX = drift.aperture(1);
               drift1.apertureY = drift.aperture(2);
            end

            calllib(trackingMethod,[libroutine 'Drift'],screen1);

        end % function TrackLibrary

    end % methods

end % classdef Screen
