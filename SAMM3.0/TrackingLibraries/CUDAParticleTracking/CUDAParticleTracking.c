#include <cuda_runtime.h>

#include "CUDAParticleTracking.h"

void SetBeam();
void GetBeam();

void TrackMarkerCUDA();
void  CopyMarkerCUDA();

void TrackDriftCUDA(DriftParameters_t drift);
void  CopyDriftCUDA(DriftParameters_t drift);

void TrackSolenoidCUDA(SolenoidParameters_t solenoid);

void TrackDipoleCUDA(DipoleParameters_t dipole);

void TrackQuadrupoleCUDA(QuadrupoleParameters_t quadrupole);

void TrackSextupoleCUDA(SextupoleParameters_t sextupole);

void TrackMultipoleCUDA(MultipoleParameters_t multipole);
void  CopyMultipoleCUDA(MultipoleParameters_t multipole);

void TrackOrbitCorrectorCUDA(OrbitCorrectorParameters_t orbitcorrector);
void  CopyOrbitCorrectorCUDA(OrbitCorrectorParameters_t orbitcorrector);

void TrackRFCavityCUDA(RFCavityParameters_t rfcavity);
void  CopyRFCavityCUDA(RFCavityParameters_t rfcavity);

void TrackBPMCUDA(BPMParameters_t bpm, double* bpmX, double* bpmY, int* nparticles);

void TrackFieldMapCUDA(FieldMapParameters_t fieldmap);

void TrackBeamlineCUDA(int n1, int n2, int npart);

BeamParameters_t* beamParams;
double* beamPS;
double* beamSpins;
double* distance;

void ParticleTracking(InitialisationParameters_t* params)
{
	cudaGetDeviceCount(&params->valid);
	params->version = "Particle tracking in CUDA version 3.1 (3 October 2013)";
}

double Precision()
{
	return (double)1e-8;
}

void InitialiseTracking(BeamParameters_t* beamParams_, double* beamPS_, double* beamSpins_, double* distance_)
{
	beamParams   = beamParams_;
	beamPS       = beamPS_;
	beamSpins    = beamSpins_;
	distance     = distance_;

	SetBeam();
}

void TrackBeamline(int n1, int n2)
{
	TrackBeamlineCUDA(n1, n2, beamParams->nparticles);
}

void DeleteBeamlineCopy()
{

}

void FinishTracking()
{
	GetBeam();
}

void TrackMarker()
{
	TrackMarkerCUDA();
}

void CopyMarker()
{
	CopyMarkerCUDA();
}

void TrackDrift(DriftParameters_t drift)
{
	TrackDriftCUDA(drift);
}

void CopyDrift(DriftParameters_t drift)
{
	CopyDriftCUDA(drift);
}

void TrackSolenoid(SolenoidParameters_t solenoid)
{
	TrackSolenoidCUDA(solenoid);
}

void TrackDipole(DipoleParameters_t dipole)
{
	TrackDipoleCUDA(dipole);
}

void TrackQuadrupole(QuadrupoleParameters_t quadrupole)
{
	TrackQuadrupoleCUDA(quadrupole);
}

void TrackSextupole(SextupoleParameters_t sextupole)
{
	TrackSextupoleCUDA(sextupole);
}

void TrackMultipole(MultipoleParameters_t multipole)
{
	TrackMultipoleCUDA(multipole);
}

void CopyMultipole(MultipoleParameters_t multipole)
{
	CopyMultipoleCUDA(multipole);
}

void TrackOrbitCorrector(OrbitCorrectorParameters_t orbitcorrector)
{
	TrackOrbitCorrectorCUDA(orbitcorrector);
}

void CopyOrbitCorrector(OrbitCorrectorParameters_t orbitcorrector)
{
	CopyOrbitCorrectorCUDA(orbitcorrector);
}

void TrackRFCavity(RFCavityParameters_t rfcavity)
{
	TrackRFCavityCUDA(rfcavity);
}

void CopyRFCavity(RFCavityParameters_t rfcavity)
{
	CopyRFCavityCUDA(rfcavity);
}

void TrackBPM(BPMParameters_t bpm, double* bpmX, double* bpmY, int* nparticles)
{
	TrackBPMCUDA(bpm, bpmX, bpmY, nparticles);
}

void TrackFieldMap(FieldMapParameters_t fieldmap)
{
	TrackFieldMapCUDA(fieldmap);
}