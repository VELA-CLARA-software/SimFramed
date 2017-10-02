% DemoTrackingLibraries1
% Tracks a bunch through a FODO lattice using different tracking libraries,
% and plots the time taken as a function of the number of particles in the
% bunch.

DefineBeamline

% Set up the particle bunch
bl.SetTrackingMethod('CParticleTracking');

[beta, tune, closedorbit] = ComputeMatchedTwiss(bl,beam);

epsilonx = 1e-6;    % horizontal emittance in metres
epsilony = 1e-10;   % vertical emittance in metres
epsilonz = 1e-4;    % longitudinal emittance in metres

nptcle   = int32(10.^(1:0.5:4));
nturns   = 1;

% Benchmark Matlab tracking
bl.SetTrackingMethod('default');
timeMatlab = zeros(size(nptcle));

for n = 1:size(nptcle,2)
    nparticles = nptcle(n); % number of particles
    beam.particles = MakeMatchedBunch(closedorbit(:,end),beta(:,:,:,end),[epsilonx, epsilony, epsilonz],nparticles);
    tStart = tic;
    beam  = bl.Track([1 nturns*length(bl.componentlist)],beam);
    timeMatlab(n) = toc(tStart);
end

figure
plot(log10(double(nptcle)),log10(timeMatlab),'-ok')
hold on
xlabel('log_{10}(number of particles)')
ylabel('log_{10}(time in seconds)')
legend('Matlab','Location','NorthWest')

% Benchmark CParticleTracking
bl.SetTrackingMethod('CParticleTracking');
timeC  = zeros(size(nptcle));

for n = 1:size(nptcle,2)
    nparticles = nptcle(n); % number of particles
    beam.particles = MakeMatchedBunch(closedorbit(:,end),beta(:,:,:,end),[epsilonx, epsilony, epsilonz],nparticles);
    tStart = tic;
    beam  = bl.Track([1 nturns*length(bl.componentlist)],beam);
    timeC(n) = toc(tStart);
end

plot(log10(double(nptcle)),log10(timeC),'-xb')
legend('Matlab','C','Location','NorthWest')

% Benchmark CUDAParticleTracking
bl.SetTrackingMethod('CUDAParticleTracking');
timeCUDA = zeros(size(nptcle));

for n = 1:size(nptcle,2)
    nparticles = nptcle(n); % number of particles
    beam.particles = MakeMatchedBunch(closedorbit(:,end),beta(:,:,:,end),[epsilonx, epsilony, epsilonz],nparticles);
    tStart = tic;
    beam = bl.Track([1 nturns*length(bl.componentlist)],beam);
    timeCUDA(n) = toc(tStart);
end

plot(log10(double(nptcle)),log10(timeCUDA),'-^m')
legend('Matlab','C','CUDA','Location','NorthWest')

bl.SetTrackingMethod('default');

set(gcf,'PaperUnits','inches')
set(gcf,'PaperPosition',[1 1 9 6])
print('-dpng','TrackingBenchmark.png')
