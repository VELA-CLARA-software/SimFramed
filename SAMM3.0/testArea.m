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
quad1.length = 0.05;
quad1.name = 'quad1';
quad1.gradient = 0.00;
dip1 = Dipole;
dip1.length = 0.4;
dip1.name = 'dip1';
dip1.angle = pi/4;
dip1.field = beam1.rigidity * dip1.angle / dip1.length;


beamlineD = Beamline;
beamlineD.componentlist = {,dip1};
beam2 = beamlineD.Track([1,1],beam1);
beam2 = beam2.particles

%[G,r] = ComputeTransferMatrix(beamlineD, [1 2], beam1, [0 0 0 0 0 0]');
