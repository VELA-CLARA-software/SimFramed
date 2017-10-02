% DemoTrackingLibraries2
% Tracks a bunch through a FODO lattice using different tracking libraries,
% and plots the emittances (which are conserved under linear, symplectic
% transport) after each turn.

% Set up the beam line
DefineBeamline

% Set up the particle bunch
bl.SetTrackingMethod('CParticleTracking');
[beta, tune, closedorbit] = ComputeMatchedTwiss(bl,beam);

epsilonx = 1e-9;    % horizontal emittance in metres
epsilony = 1e-12;   % vertical emittance in metres
epsilonz = 1e-6;    % longitudinal emittance in metres

nptcle = 1000;
ptcles = MakeMatchedBunch(closedorbit(:,end),             ...
                          beta(:,:,:,end),                ...
                          [epsilonx, epsilony, epsilonz], ...
                          nptcle);

nturns = 10;

symp   =   [ 0  1  0  0  0  0 ;...
            -1  0  0  0  0  0 ;...
             0  0  0  1  0  0 ;...
             0  0 -1  0  0  0 ;...
             0  0  0  0  0  1 ;...
             0  0  0  0 -1  0 ];

figure;

% Matlab tracking
bl.SetTrackingMethod('default');
beam.particles  = ptcles;
beam.globaltime = 0;
emittanceMatlab = zeros(3,nturns);
tic
for n = 1:nturns
    beam  = bl.Track([1 length(bl.componentlist)],beam);
    p     = beam.particles;
    p     = p - mean(p,2)*ones(1,size(p,2));
    sigma = p*p'/size(p,2);
    eps   = sort(imag(eig(sigma*symp)));
    emittanceMatlab(:,n) = eps([5 4 6]);
end
toc

subplot(3,1,1)
plot( 1e9*emittanceMatlab(1,:),'-ok')
hold on
ylabel('\epsilon_x (nm)')
legend('Matlab','Location','NorthWest')

subplot(3,1,2)
plot(1e12*emittanceMatlab(2,:),'-ok')
hold on
ylabel('\epsilon_y (pm)')

subplot(3,1,3)
plot( 1e6*emittanceMatlab(3,:),'-ok')
hold on
ylabel('\epsilon_z (\mum)')
xlabel('Turn number')

% CParticleTracking
bl.SetTrackingMethod('CParticleTracking');

beam.particles = ptcles;
beam.globaltime = 0;
emittanceC = zeros(3,nturns);
tic
for n = 1:nturns
    beam  = bl.Track([1 length(bl.componentlist)],beam);
    p     = beam.particles;
    p     = p - mean(p,2)*ones(1,size(p,2));
    sigma = p*p'/size(p,2);
    eps   = sort(imag(eig(sigma*symp)));
    emittanceC(:,n) = eps([5 4 6]);
end
toc

subplot(3,1,1)
plot( 1e9*emittanceC(1,:),'-xb')
legend('Matlab','C','Location','NorthWest')
subplot(3,1,2)
plot(1e12*emittanceC(2,:),'-xb')
subplot(3,1,3)
plot( 1e6*emittanceC(3,:),'-xb')

% CUDAParticleTracking
bl.SetTrackingMethod('CUDAParticleTracking');
beam.particles = ptcles;
beam.globaltime = 0;
emittanceCUDA = zeros(3,nturns);
tic
for n = 1:nturns
    beam  = bl.Track([1 length(bl.componentlist)],beam);
    p     = beam.particles;
    p     = p - mean(p,2)*ones(1,size(p,2));
    sigma = p*p'/size(p,2);
    eps   = sort(imag(eig(sigma*symp)));
    emittanceCUDA(:,n) = eps([5 4 6]);
end
toc

subplot(3,1,1)
plot( 1e9*emittanceCUDA(1,:),'-^m')
legend('Matlab','C','CUDA','Location','NorthWest')
subplot(3,1,2)
plot(1e12*emittanceCUDA(2,:),'-^m')
subplot(3,1,3)
plot( 1e6*emittanceCUDA(3,:),'-^m')

bl.SetTrackingMethod('default');

set(gcf,'PaperUnits','inches')
set(gcf,'PaperPosition',[1 1 5 9])
print('-dpng','TrackingPrecision.png')
