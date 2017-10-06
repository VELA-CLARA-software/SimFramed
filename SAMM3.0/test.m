clear
EBT_S01_DRIFT_01 = Drift;
EBT_S01_DRIFT_01.length=0.4;

EBT_INJ_MAG_HVCOR_01 = OrbitCorrector;
EBT_INJ_MAG_HVCOR_01.field = [0,0];
EBT_INJ_MAG_HVCOR_01.length=0.05;

EBT_S01_DRIFT_02 = Drift;
EBT_S01_DRIFT_02.length = 0.25;

EBT_INJ_DIA_WCM_01 = Drift;
EBT_INJ_DIA_WCM_01.length=0.0;

EBT_S01_DRIFT_03 = Drift;
EBT_S01_DRIFT_03.length = 0.15;

EBT_INJ_DIA_BPM_01 = BeamPositionMonitor;
EBT_INJ_DIA_BPM_01.length = 0.05;

EBT_S01_DRIFT_04 = Drift;
EBT_S01_DRIFT_04.length = 0.05;

EBT_INJ_MAG_HVCOR_02 = OrbitCorrector;
EBT_INJ_MAG_HVCOR_02.field = [0, 0];
EBT_INJ_MAG_HVCOR_02.length = 0.05;

EBT_S01_DRIFT_05 = Drift;
EBT_S01_DRIFT_05.length=0.05;

EBT_INJ_DIA_YAG_01 = Screen;

EBT_S01_DRIFT_06 = Drift;
EBT_S01_DRIFT_06.length = 0.185;

EBT_INJ_MAG_QUAD_01 = Quadrupole;
EBT_INJ_MAG_QUAD_01.length = 0.1;
EBT_INJ_MAG_QUAD_01.gradient = 0.3;

EBT_S01_DRIFT_07 = Drift;
EBT_S01_DRIFT_07.length = 0.11;

EBT_INJ_MAG_QUAD_02 = Quadrupole;
EBT_INJ_MAG_QUAD_02.length= 0.1;
EBT_INJ_MAG_QUAD_02.gradient = -0.1;

EBT_S01_DRIFT_08 = Drift;
EBT_S01_DRIFT_08.length= 0.11;

EBT_INJ_MAG_QUAD_03 = Quadrupole;
EBT_INJ_MAG_QUAD_03.length= 0.1;
EBT_INJ_MAG_QUAD_03.gradient = 0.1;

EBT_S01_DRIFT_09 = Drift;
EBT_S01_DRIFT_09.length = 0.11;

EBT_INJ_MAG_QUAD_04 = Quadrupole;
EBT_INJ_MAG_QUAD_04.length = 0.1;
EBT_INJ_MAG_QUAD_04.gradient = -0.3;

EBT_S01_DRIFT_10 = Drift;
EBT_S01_DRIFT_10.length = 0.18;

EBT_INJ_DIA_YAG_02 = Screen;

EBT_S01_DRIFT_11 = Drift;
EBT_S01_DRIFT_11.length=0.275;

EBT_INJ_MAG_HVCOR_03 = OrbitCorrector;
EBT_INJ_MAG_HVCOR_03.field=[0, 0];
EBT_INJ_MAG_HVCOR_03.length=0.05;

EBT_S01_DRIFT_12 = Drift;
EBT_S01_DRIFT_12.length=0.05;

beamline1 = Beamline;
beamline1.componentlist = {EBT_S01_DRIFT_01, ...
                          EBT_INJ_MAG_HVCOR_01, ...
                          EBT_S01_DRIFT_02, ... 
                          EBT_INJ_DIA_WCM_01, ...
                          EBT_S01_DRIFT_03, ...
                          EBT_INJ_DIA_BPM_01, ...
                          EBT_INJ_MAG_HVCOR_02, ...
                          EBT_S01_DRIFT_05, ...
                          EBT_INJ_DIA_YAG_01, ...
                          EBT_S01_DRIFT_06, ...
                          EBT_INJ_MAG_QUAD_01, ...
                          EBT_S01_DRIFT_07, ...
                          EBT_INJ_MAG_QUAD_02, ...
                          EBT_S01_DRIFT_08, ...
                          EBT_INJ_MAG_QUAD_03, ...
                          EBT_S01_DRIFT_09, ...
                          EBT_INJ_MAG_QUAD_04, ...
                          EBT_S01_DRIFT_10, ...
                          EBT_INJ_DIA_YAG_02, ...
                          EBT_S01_DRIFT_11, ...
                          EBT_INJ_MAG_HVCOR_03, ...
                          EBT_S01_DRIFT_12};
beam1 = Beam(Electron);
beam1.energy = 4.5*PhysicalUnits.MeV;
ptcle1 = [0.001,0,0,0,0,0]';
ptcle2 = [0,0,0.001,0,0,0]';
ptcle3 = [0,0,0,0,0,0]';
beam1.particles = [ptcle1, ptcle2, ptcle3];
ptcle3 = [0,0,0,0,0,0]';
DIP = Dipole;
DIP.angle = pi/4;
DIP.length=0.4;
DIP.field = beam1.rigidity * DIP.angle / DIP.length;

beamlineD = Beamline;
beamlineD.componentlist = {DIP};
beam2 = beamlineD.Track([1,1],beam1);
beam2.particles



