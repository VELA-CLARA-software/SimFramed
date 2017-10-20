#include <cuda_runtime.h>
#include <math.h>
#include <stdlib.h>

#include "CUDAParticleTracking.h"
#include "TrackCUDA.h"

extern "C" BeamParameters_t* beamParams;
extern "C" double* beamPS;
extern "C" double* beamSpins;
extern "C" double* distance;

BeamlineComponent_t* beamline;
int blocksPerGrid;
Beam* d_BeamParameters;
Particle* d_ParticleArray;

__global__ void RotateSpins(double* field, Particle* ptcle);


extern "C" __host__ int NDevices()
{
	int nval;
	int* ndevices = (int*)malloc(sizeof(int));
	cudaGetDeviceCount(ndevices);
	nval = *ndevices;
	free(ndevices);
	return nval;
}


extern "C" __host__ void SetBeam()
{
	blocksPerGrid = (beamParams->nparticles + threadsPerBlock - 1) / threadsPerBlock;


	size_t sizeptcle = beamParams->nparticles * _nptcleparams * sizeof(Particle);	
	Particle* h_ParticleArray = (Particle*)malloc(sizeptcle);

	size_t sizeBeamParameters = _nbeamparams * sizeof(Beam);
	Beam* h_BeamParameters = (Beam*)malloc(sizeBeamParameters);


	double bg = beamParams->rigidity * beamParams->charge / beamParams->mass / SpeedOfLight;


	h_BeamParameters[_globaltime] = beamParams->globaltime;
	h_BeamParameters[_mass]       = beamParams->mass;
	h_BeamParameters[_charge]     = beamParams->charge;
	h_BeamParameters[_g]          = beamParams->g;
	h_BeamParameters[_rigidity]   = beamParams->rigidity;
	h_BeamParameters[_gamma]      = sqrt(1+bg*bg);
	h_BeamParameters[_beta]       = bg/h_BeamParameters[_gamma];

	cudaMalloc((void**)&d_BeamParameters, sizeBeamParameters);
	cudaMemcpy(d_BeamParameters, h_BeamParameters, sizeBeamParameters, cudaMemcpyHostToDevice);

	free(h_BeamParameters);


	for(int n=0; n<beamParams->nparticles; n++)
	{
		h_ParticleArray[_x(n)]  = beamPS[n*6];
		h_ParticleArray[_px(n)] = beamPS[n*6 + 1];
		h_ParticleArray[_y(n)]  = beamPS[n*6 + 2];
		h_ParticleArray[_py(n)] = beamPS[n*6 + 3];
		h_ParticleArray[_ct(n)] = beamPS[n*6 + 4];
		h_ParticleArray[_dp(n)] = beamPS[n*6 + 5];
		
		if(beamParams->g != 0)
		{
			h_ParticleArray[_theta(n)] = beamSpins[n*2];
			h_ParticleArray[_phi(n)]   = beamSpins[n*2 + 1];
		}
		
		h_ParticleArray[_s(n)] = distance[n*2];
		h_ParticleArray[_f(n)] = (int)distance[n*2 + 1];
	}

	cudaMalloc((void**)&d_ParticleArray, sizeptcle);
	cudaMemcpy(d_ParticleArray, h_ParticleArray, sizeptcle, cudaMemcpyHostToDevice);
	
	free(h_ParticleArray);

}


extern "C" __host__ void GetBeam()
{
	size_t sizeptcle = beamParams->nparticles * _nptcleparams * sizeof(Particle);	
	Particle* h_ParticleArray = (Particle*)malloc(sizeptcle);

	cudaMemcpy(h_ParticleArray, d_ParticleArray, sizeptcle, cudaMemcpyDeviceToHost);
	
	for(int n=0; n<beamParams->nparticles; n++)
	{
		beamPS[n*6]     = h_ParticleArray[_x(n)];
		beamPS[n*6 + 1] = h_ParticleArray[_px(n)];
		beamPS[n*6 + 2] = h_ParticleArray[_y(n)];
		beamPS[n*6 + 3] = h_ParticleArray[_py(n)];
		beamPS[n*6 + 4] = h_ParticleArray[_ct(n)];
		beamPS[n*6 + 5] = h_ParticleArray[_dp(n)];
		
		if(beamParams->g != 0)
		{
			beamSpins[n*2]     = h_ParticleArray[_theta(n)];
			beamSpins[n*2 + 1] = h_ParticleArray[_phi(n)];
		}
		
		distance[n*2]     = h_ParticleArray[_s(n)];
		distance[n*2 + 1] = h_ParticleArray[_f(n)];
	}

	cudaFree(d_ParticleArray);
	free(h_ParticleArray);

	size_t sizeBeamParams = _nbeamparams*sizeof(Beam);
	Beam* h_BeamParameters = (Beam*)malloc(sizeBeamParams);

	cudaMemcpy(h_BeamParameters, d_BeamParameters, sizeBeamParams, cudaMemcpyDeviceToHost);

	beamParams->globaltime = h_BeamParameters[_globaltime];
	beamParams->rigidity   = h_BeamParameters[_rigidity];

	cudaFree(d_BeamParameters);
	free(h_BeamParameters);
}


extern "C" __host__ void AppendComponent(void (*trackingRoutine)(Component*, Beam*, Particle*), Component* parameterList)
{
	static BeamlineComponent_t* currentBeamlineComponent;

	BeamlineComponent_t* blcpt = (BeamlineComponent_t *)malloc(sizeof(BeamlineComponent_t));

	blcpt->trackingRoutine = trackingRoutine;
	blcpt->parameterList   = parameterList;
	blcpt->nextComponent   = NULL;

	if(beamline)
	{
		currentBeamlineComponent->nextComponent = blcpt;
		currentBeamlineComponent = blcpt;
	}
	else
	{
		beamline = blcpt;
		currentBeamlineComponent = blcpt;
	}
}


extern "C" __host__ void TrackBeamlineCUDA(int n1, int n2, int npart)
{
	BeamlineComponent_t* blcpt = beamline;

	int n = 1;

	blocksPerGrid = (npart + threadsPerBlock - 1) / threadsPerBlock;

	while(blcpt && (n<n1))
	{
		blcpt = blcpt->nextComponent;
		n++;
	}

	while(blcpt && (n<=n2))
	{
		(blcpt->trackingRoutine)<<<blocksPerGrid, threadsPerBlock>>>(blcpt->parameterList, d_BeamParameters, d_ParticleArray);
		blcpt = blcpt->nextComponent;
		n++;
	}
}
