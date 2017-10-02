
DefineBeamline

[beta tune closedorbit] = ComputeMatchedTwiss(bl,beam);

gridmax  = 0.05;
gridstep = 0.0005;

quadFmap = FieldMap;
quadFmap.length = quadF.length;
gradient = quadF.gradient;
quadFmap.gridX  = -gridmax:gridstep:gridmax;
quadFmap.gridY  = -gridmax:gridstep:gridmax;
quadFmap.gridZ  = 0:(quadF.length/20):quadF.length;
quadFmap.interpmethod = 'linear*';
quadFmap.nsteps = 200;

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

bl.SetTrackingMethod('CUDAParticleTrackingSP');
% bl.SetTrackingMethod('Matlab');

[beta1 tune1 closedorbit1] = ComputeMatchedTwiss(bl,beam);

beam.particles = MakeMatchedBunch(closedorbit(:,end),beta(:,:,:,end),[1e-8,1e-8,1e-6],100000);

tic
for n = 1:10
    beam  = bl.Track([1 length(bl.componentlist)],beam);
end
toc

unloadlibrary('CUDAParticleTrackingSP')

close all

ptcles = beam.particles;

figure
plot(ptcles(1,:),ptcles(2,:),'.b','MarkerSize',1)
xlabel('x [m]')
ylabel('p_x')

figure
plot(ptcles(3,:),ptcles(4,:),'.r','MarkerSize',1)
xlabel('y [m]')
ylabel('p_y')

figure
plot(ptcles(5,:),ptcles(6,:),'.k','MarkerSize',1)
xlabel('ct [m]')
ylabel('\delta')
