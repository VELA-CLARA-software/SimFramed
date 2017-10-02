
beam        = Beam(Positron);
beam.energy = 2.0 * PhysicalUnits.GeV;

ncells      = 16;

drift1 = Drift;
  drift1.length = 5.00;                  % metres
  
drift2 = Drift;
  drift2.length = 6.00;                  % metres
  
drift3 = Drift;
  drift3.length = 0.10;                  % metres
  
drift4 = Drift;
  drift4.length = 4.75;                  % metres
  
drift5 = Drift;
  drift5.length = 5.75;                  % metres
  
dipole1 = Dipole;
  dipole1.length = 2.50;                 % metres
  dipole1.angle  = 2*pi/ncells;          % radians
  dipole1.field  = beam.rigidity * dipole1.angle / dipole1.length; % tesla
  
quadF  = Quadrupole;
  quadF.length   = 0.25;                 % metres
  quadF.gradient = 0.45 * beam.rigidity; % tesla/metre
  quadF.aperture = [0.05, 0.05];         % metres
  
quadD  = Quadrupole;
  quadD.length   = 0.25;                 % metres
  quadD.gradient =-0.40 * beam.rigidity; % tesla/metre
  
sextF  = Sextupole;
  sextF.length   = 0.15;                 % metres
  sextF.gradient = 0.27 * beam.rigidity; % tesla/metre^2

sextD  = Sextupole;
  sextD.length   = 0.15;                 % metres
  sextD.gradient =-0.37 * beam.rigidity; % tesla/metre^2

rfcav1 = RFCavity;
  rfcav1.length   = 0.2965;              % metres
  rfcav1.voltage  = 600e3*sign(beam.species.charge); % volts
  % Note that we define a separate variable for the nominal rf frequency.
  % The actual frequency of the cavity is set below, to an harmonic of
  % the revolution frequency.
  rffreq          = 500 * PhysicalUnits.megahertz;

drift6 = Drift;
  drift6.length = drift5.length - rfcav1.length; % metres
  
bpm = cell(1,2*ncells);  
  
bl = Beamline;
  
for n = 1:ncells
    
    bl.AppendComponent(drift2);
%     bpm{2*n-1} = BeamPositionMonitor;
%     bl.AppendComponent(bpm{2*n-1});
    bl.AppendComponent(quadF);
    bl.AppendComponent(drift3);
    bl.AppendComponent(sextF);
    bl.AppendComponent(drift4);
    bl.AppendComponent(dipole1);
    bl.AppendComponent(drift1);
%     bpm{2*n}   = BeamPositionMonitor;
%     bl.AppendComponent(bpm{2*n});
    bl.AppendComponent(quadD);
    bl.AppendComponent(drift3);
    bl.AppendComponent(sextD);
    bl.AppendComponent(drift5);
    
end

bl.componentlist{end} = drift6;
bl.AppendComponent(rfcav1);

svals = bl.ComputePositions();

% Set Master Oscillator frequency before setting RF cavity harmonic and phase
MasterOscillator.SetFrequency(beam.beta*PhysicalConstants.SpeedOfLight/svals(end));
rfcav1.harmonic = floor(rffreq/MasterOscillator.GetFrequency());
rfcav1.phase    = pi - 2*pi*svals(end-1)*rfcav1.frequency/beam.beta/PhysicalConstants.SpeedOfLight;

% Check the stability of the longitudinal motion
m = ComputeTransferMatrix(bl,[1 numel(bl.componentlist)],beam,zeros(6,1));
if max(abs(eig(m(:,:,end))))-1>1e-6 % if true, then we need to reverse the polarity of the rf
    rfcav1.voltage = -1*rfcav1.voltage;
    m = ComputeTransferMatrix(bl,[1 numel(bl.componentlist)],beam,zeros(6,1));
end

% m(:,:,end)
