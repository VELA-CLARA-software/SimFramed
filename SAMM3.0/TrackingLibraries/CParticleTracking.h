
typedef struct CParticleTracking_InitialisationParameters
{
	char* version;
	int valid;
} InitialisationParameters_t;

typedef struct CParticleTracking_BeamParameters
{
	int nparticles;
	double rigidity;
	double charge;
	double mass;
	double g;
	double globaltime;
} BeamParameters_t;

typedef struct CParticleTracking_DriftParameters
{
	double length;
	double apertureX;
	double apertureY;
} DriftParameters_t;

typedef struct CParticleTracking_SolenoidParameters
{
	double length;
	double field;
	double taper;
	double apertureX;
	double apertureY;
} SolenoidParameters_t;

typedef struct CParticleTracking_DipoleParameters
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

typedef struct CParticleTracking_QuadrupoleParameters
{
	double length;
	double gradient;
	double apertureX;
	double apertureY;
} QuadrupoleParameters_t;

typedef struct CParticleTracking_SextupoleParameters
{
	double length;
	double gradient;
	double apertureX;
	double apertureY;
} SextupoleParameters_t;

typedef struct CParticleTracking_OctupoleParameters
{
	double length;
	double gradient;
	double apertureX;
	double apertureY;
} OctupoleParameters_t;

typedef struct CParticleTracking_MultipoleParameters
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

typedef struct CParticleTracking_OrbitCorrectorParameters
{
	double length;
	double fieldX;
	double fieldY;
	double apertureX;
	double apertureY;
} OrbitCorrectorParameters_t;

typedef struct CParticleTracking_RFCavityParameters
{
	double length;
	double voltage;
	double frequency;
	double phase;
	double masteroscillatorfrequency;
	double apertureX;
	double apertureY;
} RFCavityParameters_t;

typedef struct CParticleTracking_BPMParameters
{
	double apertureX;
	double apertureY;
} BPMParameters_t;

typedef struct CParticleTracking_FieldMapParameters
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
	int nsteps;
} FieldMapParameters_t;

typedef struct CParticleTracking_FieldMapGenGradientsParameters
{
	double apertureX;
	double apertureY;
	int ctablemsize;
	int ctablessize;
	int ctablelsize;
	double* ctable;
	double* mvalues;
	double* svalues;
} FieldMapGenGradientsParameters_t;

typedef struct CParticleTracking_TaylorMapParameters
{
	double apertureX;
	double apertureY;
	double length;
	int ncoeffs;
	double* coeffs;
} TaylorMapParameters_t;


__declspec(dllexport) void RotateSpins(double bx, double by, double bz, double ds, double px, double py, double dp, double* theta, double* phi);

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
__declspec(dllexport) void  CopyQuadrupole(QuadrupoleParameters_t quadrupole);

__declspec(dllexport) void TrackSextupole(SextupoleParameters_t sextupole);

__declspec(dllexport) void TrackOctupole(OctupoleParameters_t octupole);

__declspec(dllexport) void TrackMultipole(MultipoleParameters_t multipole);
__declspec(dllexport) void  CopyMultipole(MultipoleParameters_t multipole);

__declspec(dllexport) void TrackOrbitCorrector(OrbitCorrectorParameters_t orbitcorrector);
__declspec(dllexport) void  CopyOrbitCorrector(OrbitCorrectorParameters_t orbitcorrector);

__declspec(dllexport) void TrackRFCavity(RFCavityParameters_t rfcavity);
__declspec(dllexport) void  CopyRFCavity(RFCavityParameters_t rfcavity);

__declspec(dllexport) void TrackBPM(BPMParameters_t bpm, double* bpmX, double* bpmY, int* nparticles);

__declspec(dllexport) void TrackFieldMap(FieldMapParameters_t fieldmap);

__declspec(dllexport) void TrackFieldMapGenGradients(FieldMapGenGradientsParameters_t fieldmapgg);
__declspec(dllexport) void  CopyFieldMapGenGradients(FieldMapGenGradientsParameters_t fieldmapgg);

__declspec(dllexport) void TrackTaylorMap(TaylorMapParameters_t taylormap);
__declspec(dllexport) void  CopyTaylorMap(TaylorMapParameters_t taylormap);
