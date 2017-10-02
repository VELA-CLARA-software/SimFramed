classdef OrbitCorrector < handle
    % OrbitCorrector class
    % 
    % Properties:
    %   name
    %   length
    %   field
    %   aperture
    %
    % Methods:
    %   Track
    %   TrackSpin
    %   GetBField

    properties
        name      = ''; % string
        length    = 0;  % dipole length, in metres
        field     = [0 0];  % horizontal and vertical magnetic field, in tesla
        aperture  = []; % 1x2 array of elliptical aperture half-axes, in metres
    end % properties
    
    
    methods
        
        function beam = Track(orbitcorrector,beam)
            % beam2 = OrbitCorrector.Track(beam1)
            % Applies the transfer map for an orbit corrector
            % to the particles in in beam1.
            
            [x0, px0, y0, py0, ct0, dp0] = beam.GetParticles();

            beta0  = beam.beta;

            ds     = orbitcorrector.length;
            k0     = orbitcorrector.field / beam.rigidity;
            d1     = sqrt(1 + 2*dp0/beta0 + dp0.*dp0);
                        
            x1     = x0 + ds*px0./d1 - ds*ds*k0(2)./d1/2;
            px1    =         px0     - ds*k0(2);
            
            y1     = y0 + ds*py0./d1 + ds*ds*k0(1)./d1/2;
            py1    =         py0     + ds*k0(1);

            f1   = ds*(1/beta0 + dp0)./d1./d1./d1/2;
            
            c0   = ct0 + ds/beta0 - ds*(1/beta0 + dp0)./d1 - ...
                   ds*ds*f1*(k0(1)*k0(1) + k0(2)*k0(2))/3;
            
            ct1 = c0 + ...
                  ds*f1.*(k0(2)*px0 - k0(1)*py0) - ...
                  f1.*(px0.*px0 + py0.*py0);

            beam.SetParticles(x1, px1, y1, py1, ct1, dp0);

        end % function Track
        
        
        function TrackLibrary(orbitcorrector,trackingMethod,libroutine)
            
            orbitcorrector1.length    = orbitcorrector.length;
            orbitcorrector1.fieldX    = orbitcorrector.field(1);
            orbitcorrector1.fieldY    = orbitcorrector.field(2);
            
            if(~isempty(orbitcorrector.aperture))
               orbitcorrector1.apertureX = orbitcorrector.aperture(1);
               orbitcorrector1.apertureY = orbitcorrector.aperture(2);
            end
            
            calllib(trackingMethod,[libroutine 'OrbitCorrector'],orbitcorrector1);
            
        end % function TrackLibrary
        
        
        function [bx, by, bz] = GetBField(orbitcorrector,~)
            % [bx, by, bz] = OrbitCorrector.GetBField(beam)
            % Returns the magnetic field (in tesla) at the locations of
            % the particles in the beam.

            bx = orbitcorrector.field(1);
            by = orbitcorrector.field(2);
            bz = 0;
            
        end % function GetBField
        
        
        function beam = TrackSpin(orbitcorrector,beam)
            % beam2 = OrbitCorrector.Trackspin(beam1)
            % Tracks particle spins through an OrbitCorrector.

            [bx, by, bz] = orbitcorrector.GetBField(beam);
            
            Utilities.SpinRotation(beam,bx,by,bz,orbitcorrector.length);
            
        end % function TrackSpin                                                      

        
    end % methods
    
end % classdef OrbitCorrector
