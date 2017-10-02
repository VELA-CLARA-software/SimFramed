
DefineBeamline

beam.particles = [0.001 0 0 0 0 0]';

svals = bl.ComputePositions();

trajectory = zeros(6,length(svals));

trajectory(:,1) = beam.particles;

for n = 1:length(bl.componentlist)

    beam  = bl.Track([n n],beam);
    
    trajectory(:,n+1) = beam.particles;
    
end

figure
plot(svals,1e3*trajectory(1,:),'-k');
xlabel('s (m)')
ylabel('x (mm)')
