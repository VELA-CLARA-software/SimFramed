classdef Utilities
    % Utilities
    
    properties (Constant)
        
    end % properties
    
    methods (Static)
        
        function SpinRotation(beam,Bx,By,Bz,ds)
            
            [~, px, ~, py, ~, dp] = beam.GetParticles();
            
            [theta0, phi0] = beam.GetSpins();

            polsnx0      = sin(theta0).*cos(phi0);                                  % initial x component of the polarisation
            polsny0      = sin(theta0).*sin(phi0);                                  % initial y component of the polarisation
            polsnz0      = cos(theta0);                                             % initial z component of the polarisation

            species      = beam.species;
            betagammaRef = beam.momentum/species.mass/PhysicalConstants.SpeedOfLight;  % beta*gamma for a particle with the reference momentum
            gammaRef     = sqrt(betagammaRef*betagammaRef + 1);                     % relativistic gamma at the reference momentum
            betaRef      = sqrt(1 - 1/gammaRef/gammaRef);                           % relativistic beta at the reference momentum

            gamma        = (dp*betaRef + 1)*gammaRef;                               % relativistic gamma for each particle
            beta         = sqrt(1 - 1 ./gamma./gamma);                              % relativistic beta for each particle
            ptot         = beta.*gamma/betagammaRef;                                % total normalised momentum
            pz           = real(sqrt(ptot.*ptot - px.*px - py.*py));                % normalised longitudinal momentum
                                                                                    % This has to be a real number!  Something needs to be fixed...
            pdotb        = (px.*Bx + py.*By + pz.*Bz)./ptot./ptot;
            
            bParx        = pdotb.*px;
            bPary        = pdotb.*py;
            bParz        = pdotb.*pz;

            bPerpx       = Bx - bParx;
            bPerpy       = By - bPary;
            bPerpz       = Bz - bParz;

            emdg         = species.charge / species.mass ./ gamma;
            G            = (species.g - 2)/2;
            omegax       = -emdg.*((1 + G*gamma).*bPerpx + (1 + G)*bParx);
            omegay       = -emdg.*((1 + G*gamma).*bPerpy + (1 + G)*bPary);
            omegaz       = -emdg.*((1 + G*gamma).*bPerpz + (1 + G)*bParz);

            omega        = sqrt(omegax.*omegax + omegay.*omegay + omegaz.*omegaz);
            
            pdotomega    = polsnx0.*omegax + polsny0.*omegay + polsnz0.*omegaz;

            coswt        = cos(omega*ds/PhysicalConstants.SpeedOfLight);
            sinwt        = sin(omega*ds/PhysicalConstants.SpeedOfLight);

            omega(omega==0) = Inf;
            
            polsnx1      = polsnx0.*coswt + ...
                             omegax.*pdotomega.*(1-coswt)./omega./omega + ...
                             (polsnz0.*omegay - polsny0.*omegaz).*sinwt./omega;

            polsny1      = polsny0.*coswt + ...
                             omegay.*pdotomega.*(1-coswt)./omega./omega + ...
                             (polsnx0.*omegaz - polsnz0.*omegax).*sinwt./omega;

            polsnz1      = polsnz0.*coswt + ...
                             omegaz.*pdotomega.*(1-coswt)./omega./omega + ...
                             (polsny0.*omegax - polsnx0.*omegay).*sinwt./omega;

            phi1         = atan2(polsny1,polsnx1);
            theta1       = acos(polsnz1);

            beam.SetSpins(theta1, phi1);
            
        end % method SpinRotation
        
    end % methods
    
end % classdef