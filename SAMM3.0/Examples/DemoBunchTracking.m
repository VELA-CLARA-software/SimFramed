
DefineBeamline

bl.SetTrackingMethod('CParticleTracking');

[beta, tune, closedorbit] = ComputeMatchedTwiss(bl,beam);

emittance  = [0.34 0.03 32.0] * PhysicalUnits.micrometre;
nparticles = 1000; % number of particles

beam.particles = MakeMatchedBunch(closedorbit(:,end),beta(:,:,:,end),emittance,nparticles);

beam.particles(:,1) = [sqrt(beta(1,1,1,1)*emittance(1)) 0 0 0 0 0 ]';
beam.particles(:,2) = [0 0 sqrt(beta(3,3,2,1)*emittance(2)) 0 0 0 ]';
beam.particles(:,3) = [0 0 0 0 sqrt(beta(5,5,3,1)*emittance(3)) 0 ]';
                       
beam.spins(1,:) = 0.1*ones(1,nparticles);
beam.spins(2,:) = (pi/2)*ones(1,nparticles);

close all
figure(1)
subplot(2,3,1); hold on; box on; axes(1) = gca;
plot(1000*beam.particles(1,:),1000*beam.particles(2,:),'.k','MarkerSize',3)
set(gca,'xlim',[-3000 3000]*sqrt(beta(1,1,1,1)*emittance(1)),'ylim',[-3000 3000]*sqrt(beta(2,2,1,1)*emittance(1)))
xlabel('x (mm)'); ylabel('p_x (10^{-3})'); title('initial phase space')

subplot(2,3,2); hold on; box on; axes(2) = gca;
plot(1000*beam.particles(3,:),1000*beam.particles(4,:),'.k','MarkerSize',3)
set(gca,'xlim',[-3000 3000]*sqrt(beta(3,3,2,1)*emittance(2)),'ylim',[-3000 3000]*sqrt(beta(4,4,2,1)*emittance(2)))
xlabel('y (mm)'); ylabel('p_y (10^{-3})'); title('initial phase space')

subplot(2,3,3); hold on; box on; axes(3) = gca;
plot(1000*beam.particles(5,:),1000*beam.particles(6,:),'.k','MarkerSize',3)
set(gca,'xlim',[-3000 3000]*sqrt(beta(5,5,3,1)*emittance(3)),'ylim',[-3000 3000]*sqrt(beta(6,6,3,1)*emittance(3)))
xlabel('z (mm)'); ylabel('\delta (10^{-3})'); title('initial phase space')

nturns = round(1/tune(end,3)); % number of turns
orbit  = zeros(6,3,nturns);
spins  = zeros(2,nparticles,nturns);

tic
for n = 1:nturns

    beam  = bl.Track([1 length(bl.componentlist)],beam);
    
    orbit(:,:,n) = beam.particles(:,1:3);
    spins(:,:,n) = beam.spins(:,:);
    
end
toc

bl.SetTrackingMethod('default');

figure(1)
for n = 1:3
    
    subplot(2,3,n)
    plot(1000*squeeze(orbit(2*n-1,n,:)),1000*squeeze(orbit(2*n,n,:)),'.r','MarkerSize',5)
    
    subplot(2,3,n+3); hold on; box on
    plot(1000*beam.particles(2*n-1,:),1000*beam.particles(2*n,:),'.k','MarkerSize',3)
    set(gca,'xlim',get(axes(n),'xlim'),'ylim',get(axes(n),'ylim'))
    xlabel(get(get(axes(n),'xlabel'),'string')); ylabel(get(get(axes(n),'ylabel'),'string')); title('final phase space')
    
end

% set(gcf,'PaperUnits','inches')
% set(gcf,'PaperPosition',[1 1 9 5])
% print('-dpng','BunchTrackingOrbit.png','-r600')

sx = mean(squeeze(sin(spins(1,:,:)).*cos(spins(2,:,:))));
sy = mean(squeeze(sin(spins(1,:,:)).*sin(spins(2,:,:))));
sz = mean(squeeze(cos(spins(1,:,:))));

figure(2)
subplot(2,1,1); hold on; box on;
plot(1:nturns, sx, '-r')
plot(1:nturns, sy, '-b')
plot(1:nturns, sz, '-k')
legend('\langleS_x\rangle','\langleS_y\rangle','\langleS_z\rangle')
xlabel('turn number'); ylabel('\langleS_i\rangle')

subplot(2,1,2)
plot(1:nturns, sqrt(sx.*sx + sy.*sy + sz.*sz))
xlabel('turn number'); ylabel('\Sigma\langleS_i\rangle^2')

% set(gcf,'PaperUnits','inches')
% set(gcf,'PaperPosition',[1 1 9 5])
% print('-dpng','BunchTrackingSpins.png','-r600')
