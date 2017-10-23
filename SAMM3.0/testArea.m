beam1 = Beam(Electron);
beam1.energy = 4.5*PhysicalUnits.MeV;
ptcle1 = [0.001,0,0,0,0,0]';
ptcle2 = [0,0,0.001,0,0,0]';
ptcle3 = [0,0,0,0,0,0]';
beam1.particles = [ptcle1, ptcle2, ptcle3];

drift = Drift;
drift.length = 0.1;
drift.name = 'drift1';
quad1 = Quadrupole;
quad1.length = 0.0;
quad1.name = 'quad1';
quad1.gradient = 0.05;


beamlineD = Beamline;
beamlineD.componentlist = {drift,quad1};
%beam2 = beamlineD.Track([1,2],beam1);
%beam2 = beam2.particles;

[G,r] = ComputeTransferMatrix(beamlineD, [1 2], beam1, [0 0 0 0 0 0]');
