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
    
    properties
        name     = ''; % string
        length   = 0;  % in metres
        aperture = []; % 1x2 array of elliptical aperture half-axes, in metres
    end % properties
    
    methods

        function beam = Track(screen,beam)
            % beam2 = Drift.Track(beam1)
            % Applies the transfer map for a drift to the particles
            % in beam1.
            
            [x0, px0, y0, py0, ct0, dp0] = beam.GetParticles();
             
            beam.SetParticles(x0, px0, y0, py0, ct0, dp0);
                        
        end % function Track
        
        
        function TrackLibrary(drift,trackingMethod,libroutine)
            
            drift1.length = drift.length;
            if(~isempty(drift.aperture))
               drift1.apertureX = drift.aperture(1);
               drift1.apertureY = drift.aperture(2);
            end
            
            calllib(trackingMethod,[libroutine 'Drift'],drift1);
            
        end % function TrackLibrary
                
    end % methods
    
end % classdef Screen