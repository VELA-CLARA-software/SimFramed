classdef Octupole < handle
    % Octupole class
    % 
    % Properties:
    %   name
    %   length
    %   gradient
    %   aperture
    %
    % Methods:
    %   Track
    %   TrackSpin
    %   GetBField
    
    
    properties
        length   = 0;  % in metres
        gradient = 0;  % in tesla/metre^2
        aperture = []; % 1x2 array of elliptical aperture half-axes, in metres
    end % properties

    methods
        
        function beam = Track(octupole,beam)
            % beam2 = Octupole.Track(beam1)
            % Applies the transfer map for a octupole to the particles
            % in beam1.
                        
            [x0, px0, y0, py0, ct0, dp0] = beam.GetParticles();
            
            beta0  = beam.beta;
            
            ds = octupole.length;
            k3 = octupole.gradient / beam.rigidity; % normalised gradient

            % First apply a drift through ds/2
            d1  = sqrt(1 - px0.*px0 - py0.*py0 + 2*dp0/beta0 + dp0.*dp0);

            x1  = x0  + ds*px0./d1/2;
            y1  = y0  + ds*py0./d1/2;
            ct1 = ct0 + ds*(1 - (1 + beta0*dp0)./d1)/beta0/2;
            
            % Next, apply a octupole 'kick'
            px1 = px0 - (x1.*x1.*x1 - 3*x1.*y1.*y1)*k3*ds/6;
            py1 = py0 + (3*x1.*x1.*y1 - y1.*y1.*y1)*k3*ds/6;
            
            % Finally, apply a second drift through ds/2
            d1  = sqrt(1 - px1.*px1 - py1.*py1 + 2*dp0/beta0 + dp0.*dp0);

            x2  = x1  + ds*px1./d1/2;
            y2  = y1  + ds*py1./d1/2;
            ct2 = ct1 + ds*(1 - (1 + beta0*dp0)./d1)/beta0/2;
            
            beam.SetParticles(x2, px1, y2, py1, ct2, dp0);
            
        end % function Track
        
        
        function TrackLibrary(octupole,trackingMethod,libroutine)
            
            octupole1.length   = octupole.length;
            octupole1.gradient = octupole.gradient;
            if(~isempty(octupole.aperture))
               octupole1.apertureX = octupole.aperture(1);
               octupole1.apertureY = octupole.aperture(2);
            end
            
            calllib(trackingMethod,[libroutine 'Octupole'],octupole1);
            
        end % function TrackLibrary
        
        
        function [bx, by, bz] = GetBField(octupole,beam)
            % [bx, by, bz] = Octupole.GetBField(beam)
            % Returns the magnetic field (in tesla) at the locations of
            % the particles in the beam.
            
            [x, ~, y] = beam.GetParticles();
            
            bx = octupole.gradient * (x.*x.*x - 3*x.*y.*y) / 6;
            by = octupole.gradient * (3*x.*x.*y - y.*y.*y) / 6;
            bz = 0;
            
        end % function GetBField
        
        
        function beam = TrackSpin(octupole,beam)
            % beam2 = Octupole.Trackspin(beam1)
            % Tracks particle spins through an Octupole.
            
            [bx, by, bz] = octupole.GetBField(beam);
            
            Utilities.SpinRotation(beam,bx,by,bz,octupole.length);

        end % function TrackSpin

        
    end % methods
    
end % classdef Octupole
