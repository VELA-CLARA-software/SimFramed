
DefineBeamline

f0 = rfcav1.frequency;

% bl.SetTrackingMethod('CParticleTracking');

I1 = ComputeRadiationIntegrals(bl,beam);
alphap = I1/svals(end);
etap   = alphap - 1/beam.gamma^2;

dpmax    = 0.006;
npts     = 5;
dpvals   = zeros(1,npts);
tunevals = zeros(2,npts);

for n = 1:npts

    dp = dpmax * (2*(n-1)/(npts-1) - 1);
    
    rfcav1.frequency = f0*(1 - etap*dp);

    [~, tune, closedorbit] = ComputeMatchedTwiss(bl,beam);

    dpvals(n) = closedorbit(6,1);
    tunevals(:,n) = tune(end,1:2)';
    
end

bl.SetTrackingMethod('default');

qx = polyfit(dpvals,tunevals(1,:),1);
qy = polyfit(dpvals,tunevals(2,:),1);

figure
subplot(2,1,1)
plot(dpvals,tunevals(1,:),'.-k')
axis tight
ylabel('\nu_x')
title(['\nu_x = ',num2str(qx(end),'%6.4f'),'; \xi_x = ',num2str(qx(end-1),'%5.3f')])
subplot(2,1,2)
plot(dpvals,tunevals(2,:),'.-r')
axis tight
xlabel('\delta')
ylabel('\nu_y')
title(['\nu_y = ',num2str(qy(end),'%6.4f'),'; \xi_y = ',num2str(qy(end-1),'%5.3f')])

set(gcf,'PaperUnits','inches')
set(gcf,'PaperPosition',[1 1 9 5])
print('-dpng','Chromaticity.png','-r600')
