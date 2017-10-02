classdef Drift < handle
    % Drift class
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

        function beam = Track(drift,beam)
            % beam2 = Drift.Track(beam1)
            % Applies the transfer map for a drift to the particles
            % in beam1.
            
            [x0, px0, y0, py0, ct0, dp0] = beam.GetParticles();
            
            beta0  = beam.beta;

            ds  = drift.length;
            
            d1  = sqrt(1 - px0.*px0 - py0.*py0 + 2*dp0/beta0 + dp0.*dp0);

            x1  = x0  + ds*px0./d1;
            y1  = y0  + ds*py0./d1;
            ct1 = ct0 + ds*(1 - (1 + beta0*dp0)./d1)/beta0;
            
            beam.SetParticles(x1, px0, y1, py0, ct1, dp0);
                        
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
    
end % classdef Drift
