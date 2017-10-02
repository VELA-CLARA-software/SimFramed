classdef Beam < handle
    % Beam class
    % 
    % Properties:
    %   species
    %   bunchcharge
    %   distance
    %   globaltime
    %   particles
    %   spins
    %   energy
    %   momentum
    %   rigidity
    %   beta ('get' only)
    %   gamma ('get' only)
    %
    % Methods:
    %   SetParticles
    %   GetParticles
    %   SetSpins
    %   GetSpins
    
    properties
        species;             % particle type, e.g. electron, positron, proton
        bunchcharge = 0;     % bunch charge in coulombs
        distance    = [];    % 2xN array, recording total distance travelled by each particle
        globaltime  = 0;     % global time in seconds
        spins       = [];    % 2xN array, polar angles of spin vector of N particles
    end % properties
    
    properties (Access=private)
        phasespace;          % phase space coordinates of particles
        brho;                % beam rigidity in tesla metres
    end % properties (private)

    properties (Dependent=true)
        particles;           % 6xN array, phase space variables of N particles
        energy;              % reference energy in joules
        momentum;            % reference momentum in kilogram metres/second
        rigidity;            % beam rigidity in tesla metres
        beta;                % relativistic beta (v/c)
        gamma;               % relativistic gamma factor
    end % properties (dependent)

    
    methods

        function beam = Beam(species)
            % initialise a beam with specified particle species
            beam.species = species;
        end % function Beam
        

        function set.energy(beam,energy)
            spec = beam.species;
            beam.rigidity = sqrt(energy^2/PhysicalConstants.SpeedOfLight^2 ...
                            - spec.mass^2*PhysicalConstants.SpeedOfLight^2)/spec.charge;
        end % energy set method
        
        function energy = get.energy(beam)
            spec = beam.species;
            energy = sqrt((beam.brho * spec.charge)^2*PhysicalConstants.SpeedOfLight^2 ...
                                        + spec.mass^2*PhysicalConstants.SpeedOfLight^4);
        end % energy get method
        
        
        function set.momentum(beam,momentum)
            spec = beam.species;
            beam.rigidity = momentum/spec.charge;
        end % momentum set method
        
        function momentum = get.momentum(beam)
            spec = beam.species;
            momentum = beam.brho * spec.charge;
        end % momentum get method

        
        function set.rigidity(beam,rigidity)
            
            spec = beam.species;
            rigidity = abs(rigidity)*sign(spec.charge);
            
            if ~isempty(beam.phasespace) && beam.brho~=0
                
                P0 = beam.brho * spec.charge;
                P1 = rigidity  * spec.charge;
                
                E0 = sqrt(P0^2*PhysicalConstants.SpeedOfLight^2 ...
                 + spec.mass^2*PhysicalConstants.SpeedOfLight^4);
                E1 = sqrt(P1^2*PhysicalConstants.SpeedOfLight^2 ...
                 + spec.mass^2*PhysicalConstants.SpeedOfLight^4);
             
                b0 = P0*PhysicalConstants.SpeedOfLight/E0;
                b1 = P1*PhysicalConstants.SpeedOfLight/E1;
                
                beam.phasespace([2,4],:) = beam.phasespace([2,4],:)*P0/P1;
                
                beam.phasespace(6,:) = (beam.phasespace(6,:) + 1/b0)*P0/P1 - 1/b1;
                
            end
            
            beam.brho = rigidity;
            
        end % rigidity set method
        
        function rigidity = get.rigidity(beam)
            rigidity = beam.brho;
        end % rigidity get method
        
        
        function beta = get.beta(beam)
            spec = beam.species;
            bg = beam.brho * spec.charge / spec.mass / PhysicalConstants.SpeedOfLight;
            beta = bg/sqrt(1+bg*bg);
        end % beta get method
        
        
        function gamma = get.gamma(beam)
            spec = beam.species;
            bg = beam.brho * spec.charge / spec.mass / PhysicalConstants.SpeedOfLight;
            gamma = sqrt(1+bg*bg);
        end % gamma get method
        
        
        function set.particles(beam,particles)
            beam.phasespace = particles;
            beam.distance   = [ zeros(1,size(particles,2));...
                                ones( 1,size(particles,2)) ];
        end % particles set method
        
        function ps = get.particles(beam)
            ps = beam.phasespace;
        end % particles get method
        
        
        function val = SetParticles(beam,x,px,y,py,ct,dp)
            % Beam.SetParticles(x,px,y,py,ct,dp)

            beam.phasespace(:,beam.distance(2,:)~=0) = [x; px; y; py; ct; dp];
            val = 0; 
        end % SetParticles method
        
        function [x, px, y, py, ct, dp] = GetParticles(beam)
            % [x, px, y, py, ct, dp] = Beam.GetParticles()
            
            inbeam = find(beam.distance(2,:));
            x      = beam.phasespace(1,inbeam);
            px     = beam.phasespace(2,inbeam);
            y      = beam.phasespace(3,inbeam);
            py     = beam.phasespace(4,inbeam);
            ct     = beam.phasespace(5,inbeam);
            dp     = beam.phasespace(6,inbeam);
        end % GetParticles method
        
                
        function val = SetSpins(beam,theta,phi)
            % Beam.SetSpins(theta,phi)

            beam.spins(:,beam.distance(2,:)~=0) = [theta; phi];
            val = 0;
        end % SetSpins method
                
        function [theta, phi] = GetSpins(beam)
            % [theta, phi] = Beam.GetSpins()
            
            inbeam = find(beam.distance(2,:)~=0);            
            theta  = beam.spins(1,inbeam);
            phi    = beam.spins(2,inbeam);
        end % GetSpins method
        
    end % methods


end % classdef Beam
