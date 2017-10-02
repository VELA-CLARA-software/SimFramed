
DefineBeamline

% bl.SetTrackingMethod('CParticleTracking')

f0 = MasterOscillator.GetFrequency();
MasterOscillator.SetFrequency(f0*0.9999);

[beta, tune, closedorbit] = ComputeMatchedTwiss(bl,beam);
[I1, I2, I3, I4, I5] = ComputeRadiationIntegrals(bl,beam);

% bl.SetTrackingMethod('default')

u0 = beam.species.Cgamma*(beam.energy^4)*I2/2/pi;
t0 = svals(end)/(beam.beta*PhysicalConstants.SpeedOfLight);
tau = 2*t0*beam.energy/u0;
jx = 1 - I4/I2;
jz = 2 + I4/I2;
gamma = beam.gamma;
sigd = gamma*sqrt(beam.species.Cq * I3/jz/I2);
sigz = I1*sigd/2/pi/tune(end,3);
eps0 = gamma*gamma*beam.species.Cq * I5/jx/I2;

fid = fopen('RingAnalysis.dat','w');

fprintf(fid,'Beam energy                = %4.6g GeV\n',beam.energy/PhysicalUnits.GeV);
fprintf(fid,'Circumference              = %4.6g m\n',svals(end));
fprintf(fid,'Revolution frequency       = %4.6g kHz\n',1/t0/PhysicalUnits.kilohertz);
fprintf(fid,'Horizontal tune            = %4.6g\n',tune(end,1));
fprintf(fid,'Vertical tune              = %4.6g\n',tune(end,2));
fprintf(fid,'Longitudinal tune          = %4.6g\n',tune(end,3));
fprintf(fid,'RF frequency               = %4.6g MHz\n',rfcav1.frequency/PhysicalUnits.megahertz);
fprintf(fid,'Harmonic number            = %4.6g\n',rfcav1.frequency/(PhysicalConstants.SpeedOfLight/svals(end)));
fprintf(fid,'Momentum compaction factor = %4.3g\n',I1/svals(end));
fprintf(fid,'Energy loss per turn       = %4.3g keV\n',u0/PhysicalUnits.keV);
fprintf(fid,'Horizontal damping partition number = %4.3g\n',jx);
fprintf(fid,'Horizontal damping time    = %4.3g ms\n',tau/jx/PhysicalUnits.millisecond);
fprintf(fid,'Vertical damping time      = %4.3g ms\n',tau/PhysicalUnits.millisecond);
fprintf(fid,'Longitudinal damping time  = %4.3g ms\n',tau/(3-jx)/PhysicalUnits.millisecond);
fprintf(fid,'Natural rms energy spread  = %4.3g\n',sigd);
fprintf(fid,'Natural rms bunch length   = %4.3g mm\n',1000*sigz);
fprintf(fid,'Natural emittance          = %4.3g nm\n',eps0/PhysicalUnits.nanometre);

fclose(fid);

beta1 = permute(beta,[4 1 2 3]);

close all
figure

subplot(3,1,1)
plot(svals, 1000*closedorbit(1,:), '-k');
hold on
plot(svals, 1000*closedorbit(3,:), '-r');
axis([0 svals(end) -inf inf])
xlabel('s [m]');
ylabel('closed orbit [mm]');

subplot(3,1,2)
plot(svals, beta1(:,1,1,1), '-k');
hold on
plot(svals, beta1(:,3,3,2), '-r');
axis([0 svals(end) -inf inf])
xlabel('s [m]');
ylabel('\beta [m]');

subplot(3,1,3)
plot(svals, beta1(:,1,6,3)./beta1(:,6,6,3), '-k');
axis([0 svals(end) -inf inf])
xlabel('s [m]');
ylabel('\eta [m]');

set(gcf,'PaperUnits','inches')
set(gcf,'PaperPosition',[1 1 9 5])
print('-dpng','RingAnalysis.png','-r600')
