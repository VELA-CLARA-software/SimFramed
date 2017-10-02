
DefineBeamline

[beta tune closedorbit] = ComputeMatchedTwiss(bl,beam);

gridmax  = 0.05;
gridstep = 0.001;

quadFmap = FieldMap;
quadFmap.length = quadF.length;
gradient = quadF.gradient;
quadFmap.gridX  = -gridmax:gridstep:gridmax;
quadFmap.gridY  = -gridmax:gridstep:gridmax;
quadFmap.gridZ  = 0:(quadF.length/20):quadF.length;
quadFmap.interpmethod = 'linear*';
quadFmap.nsteps = 100;

for nx = 1:size(quadFmap.gridX,2)
   
    x = -gridmax + gridstep*(nx-1);
    
    for ny = 1:size(quadFmap.gridY,2)

        y = -gridmax + gridstep*(ny-1);
        
        for nz = 1:size(quadFmap.gridZ,2)

            z = (quadF.length/20)*(nz-1);
            
            quadFmap.Bx(nx,ny,nz) = quadF.gradient * y;
            quadFmap.By(nx,ny,nz) = quadF.gradient * x;
            quadFmap.Bz(nx,ny,nz) = 0;
            
        end
        
    end
    
end

bl.componentlist{4} = quadFmap;

ptcle1 = closedorbit(:,1) + [0.001 0.001 0 0 0 0]';
ptcle2 = closedorbit(:,1) + [0 0 0.001 0.001 0 0]';
ptcle3 = closedorbit(:,1) + [0 0 0 0 0 0.001]';
beam.particles = [ptcle1 ptcle2 ptcle3];

nturns = 100;
x  = zeros(1,nturns);
px = zeros(1,nturns);
y  = zeros(1,nturns);
py = zeros(1,nturns);
ct = zeros(1,nturns);
dp = zeros(1,nturns);

bl.SetTrackingMethod('CUDAParticleTracking');

tic
for n = 1:nturns
    beam  = bl.Track([1 length(bl.componentlist)],beam);
    x(n)  = beam.particles(1,1);
    px(n) = beam.particles(2,1);
    y(n)  = beam.particles(3,2);
    py(n) = beam.particles(4,2);
    ct(n) = beam.particles(5,3);
    dp(n) = beam.particles(6,3);
end
toc

unloadlibrary('CUDAParticleTracking')

if nturns>10

    close all

    figure
    plot(x,px,'.b')
    xlabel('x [m]')
    ylabel('p_x')

    figure
    plot(y,py,'.r')
    xlabel('y [m]')
    ylabel('p_y')

    figure
    plot(ct,dp,'.k')
    xlabel('ct [m]')
    ylabel('\delta')

end
