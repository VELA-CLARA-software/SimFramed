classdef RFAcceleratingStructure < handle
    % RFAcceleratingStructure class
    %
    % RFAcceleratingStructure:
    %   name
    %   structuretype
    %   length
    %   ncell
    %   voltage
    %   harmonic
    %   phase
    %   globalclock
    %   aperture
    %
    % Methods:
    %   Track

    properties
        name          = ''; % string
        structuretype = 'TravellingWave'; % 'TravellingWave' or 'StandingWave'
        length        = 0;  % in metres
        ncell         = 1;  % number of cells
        voltage       = 0;  % in volts
        harmonic      = 1;  % frequency (in units of master oscillator frequency)
        phase         = pi; % relative to master oscillator, in radians
        globalclock   = 1;  % 1 = synchronise with global clock
        aperture      = []; % 1x2 array of elliptical aperture half-axes, in metres
    end % properties

    properties (Dependent=true)
        frequency;
    end % properties (dependent)

    methods

        function set.frequency(rfacceleratingstructure,f)
            f1 = rfacceleratingstructure.harmonic * MasterOscillator.GetFrequency();
            MasterOscillator.SetFrequency(MasterOscillator.GetFrequency()*f/f1);
        end

        function f = get.frequency(rfacceleratingstructure)
            f = rfacceleratingstructure.harmonic * MasterOscillator.GetFrequency();
        end

        function beam = Track(rfstructure,beam)
            % beam2 = RFAcceleratingStructure.Track(beam1)
            % Applies the dynamical map for a linac structure to the particles
            % in beam1.  The reference energy is changed by qV*cos(phase).

            [x0, px0, y0, py0, ct0, dp0] = beam.GetParticles();

            mofreq = MasterOscillator.GetFrequency();

            nc     = rfstructure.ncell;
            f      = rfstructure.harmonic*mofreq;
            L      = rfstructure.length / nc;
            cosphi = cos(rfstructure.phase);

            beam.globaltime = beam.globaltime - floor(beam.globaltime*mofreq)/mofreq;

            gt = beam.globaltime * rfstructure.globalclock;

            dE = beam.species.charge * rfstructure.voltage * cosphi / nc;
            E0 = beam.energy;
            E1 = E0 + dE;

            if abs(dE/E0)>1e-6

                if strcmp(rfstructure.structuretype,'TravellingWave')

                    logE = log(1 + dE/E0);

                    r11 = 1 - logE/2;
                    r12 = L * logE * E0/dE;
                    r21 = -dE * logE / 4 / L / E1;
                    r22 = E0 * (2 + logE) / 2 / E1;

                elseif strcmp(rfstructure.structuretype,'StandingWave')

                    L1 = PhysicalConstants.SpeedOfLight / 2 / f;
                    if abs(L/L1 - 1)>1e-2
                        error('RFAcceleratingStructure:BadLength',...
                            'RFAcceleratingStructure.length should be c/2f.')
                    end

                    a = log(E1/E0)/sqrt(8)/cosphi;
                    Eprime = (E1 - E0)/L;

                    r11 = cos(a) - sqrt(2)*cosphi*sin(a);
                    r12 = sqrt(8)*cosphi*sin(a)*E0/Eprime;
                    r21 =-(cosphi/sqrt(2) + 1/sqrt(8)/cosphi)*sin(a)*Eprime/E1;
                    r22 = (cos(a) + sqrt(2)*cosphi*sin(a))*E0/E1;

                else

                    error('RFAcceleratingStructure:UnrecognisedType',...
                        'RFAcceleratingStructure.structuretype should be StandingWave or TravellingWave.');

                end

            else

                r11 = 1;
                r12 = L * (1 - dE/2/E0);
                r21 = 0;
                r22 = 1 - dE/E0;

            end

            for n = 1:nc

                % First, apply a drift map through L/2
                % to the longitudinal coordinate
                beta0  = beam.beta;
                d1     = sqrt(1 - px0.*px0 - py0.*py0 + 2*dp0/beta0 + dp0.*dp0);
                ct0    = ct0 + L*(1 - (1 + beta0*dp0)./d1)/beta0/2;

                % Now apply the RF structure map to the transverse variables
                % and the momentum deviation
                x1  = r11*x0 + r12*px0;
                px0 = r21*x0 + r22*px0;
                x0  = x1;

                y1  = r11*y0 + r12*py0;
                py0 = r21*y0 + r22*py0;
                y0  = y1;

                P0  = beam.momentum;
                Edc = (dp0 + 1/beta0)*P0;

                beam.energy = E1;
                beta0 = beam.beta;
                P0    = beam.momentum;

                t = gt - ct0/PhysicalConstants.SpeedOfLight;
                Edc = Edc + ...
                    beam.species.charge * rfstructure.voltage * cos(2*pi*f*t + rfstructure.phase) / nc / PhysicalConstants.SpeedOfLight;

                dp0 = Edc/P0 - 1/beta0;

                % Finally, apply a drift map through L/2
                % to the longitudinal coordinate
                d1     = sqrt(1 - px0.*px0 - py0.*py0 + 2*dp0/beta0 + dp0.*dp0);
                ct0    = ct0 + L*(1 - (1 + beta0*dp0)./d1)/beta0/2;

            end

            % Set the new values for the dynamical variables
            beam.SetParticles(x0, px0, y0, py0, ct0, dp0);

        end % function Track

        function TrackLibrary(rfstructure,trackingMethod,libroutine)

            rfstructure1.length   = rfstructure.length;
            rfstructure1.voltage  = rfstructure.voltage;
            rfstructure1.harmonic = rfstructure.harmonic;
            rfstructure1.phase    = rfstructure.phase;

            rfstructure1.masteroscillatorfrequency = MasterOscillator.GetFrequency();

            if(~isempty(rfstructure.aperture))
               rfstructure1.apertureX = rfstructure.aperture(1);
               rfstructure1.apertureY = rfstructure.aperture(2);
            end

            calllib(trackingMethod,[libroutine 'RFAcceleratingStructure'],rfstructure1);

        end % function TrackLibrary


    end % methods

end % classdef RFAcceleratingStructure
