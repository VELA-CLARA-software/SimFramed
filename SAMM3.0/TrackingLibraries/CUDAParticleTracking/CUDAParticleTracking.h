#define SpeedOfLight (double)2.99792458e8

typedef struct CUDAParticleTracking_InitialisationParameters
{
	char* version;
	int valid;
} InitialisationParameters_t;

typedef struct CUDAParticleTracking_BeamParameters
{
	int nparticles;
	double rigidity;
	double charge;
	double mass;
	double g;
	double globaltime;
} BeamParameters_t;

typedef struct CUDAParticleTracking_DriftParameters
{
	double length;
	double apertureX;
	double apertureY;
} DriftParameters_t;

typedef struct CUDAParticleTracking_SolenoidParameters
{
	double length;
	double field;
	double taper;
	double apertureX;
	double apertureY;
} SolenoidParameters_t;

typedef struct CUDAParticleTracking_DipoleParameters
{
	double length;
	double field;
	double curvature;
	double gradient;
	double hgap;
	double e1;
	double fint1;
	double e2;
	double fint2;
	double apertureX;
	double apertureY;
} DipoleParameters_t;

typedef struct CUDAParticleTracking_QuadrupoleParameters
{
	double length;
	double gradient;
	double apertureX;
	double apertureY;
} QuadrupoleParameters_t;

typedef struct CUDAParticleTracking_SextupoleParameters
{
	double length;
	double gradient;
	double apertureX;
	double apertureY;
} SextupoleParameters_t;

typedef struct CUDAParticleTracking_MultipoleParameters
{
	double angle;
	double curvature;
	int ncpts;
	double apertureX;
	double apertureY;
	double* fieldindex;
	double* bnL;
	double* anL;
} MultipoleParameters_t;

typedef struct CUDAParticleTracking_OrbitCorrectorParameters
{
	double length;
	double fieldX;
	double fieldY;
	double apertureX;
	double apertureY;
} OrbitCorrectorParameters_t;

typedef struct CUDAParticleTracking_RFCavityParameters
{
	double length;
	double voltage;
	double frequency;
	double phase;
	double apertureX;
	double apertureY;
	double masteroscillatorfrequency;
} RFCavityParameters_t;

typedef struct CUDAParticleTracking_BPMParameters
{
	double apertureX;
	double apertureY;
} BPMParameters_t;

typedef struct CUDAParticleTracking_FieldMapParameters
{
	double length;
	double apertureX;
	double apertureY;
	int gridXsize;
	int gridYsize;
	int gridZsize;
	double* gridX;
	double* gridY;
	double* gridZ;
	double* Bx;
	double* By;
	double* Bz;
	int nsteptcle;
} FieldMapParameters_t;


__declspec(dllexport) void ParticleTracking(InitialisationParameters_t* params);

__declspec(dllexport) double Precision();

__declspec(dllexport) void InitialiseTracking(BeamParameters_t* beamParams_, double* beamPS_, double* beamSpins_, double* distance_);

__declspec(dllexport) void TrackBeamline(int n1, int n2);

__declspec(dllexport) void DeleteBeamlineCopy();

__declspec(dllexport) void FinishTracking();

__declspec(dllexport) void TrackMarker();
__declspec(dllexport) void  CopyMarker();

__declspec(dllexport) void TrackDrift(DriftParameters_t drift);
__declspec(dllexport) void  CopyDrift(DriftParameters_t drift);

__declspec(dllexport) void TrackSolenoid(SolenoidParameters_t solenoid);

__declspec(dllexport) void TrackDipole(DipoleParameters_t dipole);

__declspec(dllexport) void TrackQuadrupole(QuadrupoleParameters_t quadrupole);

__declspec(dllexport) void TrackSextupole(SextupoleParameters_t sextupole);

__declspec(dllexport) void TrackMultipole(MultipoleParameters_t multipole);
__declspec(dllexport) void  CopyMultipole(MultipoleParameters_t multipole);

__declspec(dllexport) void TrackOrbitCorrector(OrbitCorrectorParameters_t orbitcorrector);
__declspec(dllexport) void  CopyOrbitCorrector(OrbitCorrectorParameters_t orbitcorrector);

__declspec(dllexport) void TrackRFCavity(RFCavityParameters_t rfcavity);
__declspec(dllexport) void  CopyRFCavity(RFCavityParameters_t rfcavity);

__declspec(dllexport) void TrackBPM(BPMParameters_t bpm, double* bpmX, double* bpmY, int* nparticles);

__declspec(dllexport) void TrackFieldMap(FieldMapParameters_t fieldmap);
