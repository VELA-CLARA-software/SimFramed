#include <cuda_runtime.h>
#include <math.h>

#include "CUDAParticleTracking.h"
#include "TrackCUDA.h"

#define _nparams 3
#define _length 0
#define _apertureX 1
#define _apertureY 2

extern Beam* d_BeamParameters;
extern Particle* d_ParticleArray;
extern int blocksPerGrid;

__global__ void DriftMap(Component* drift, Beam* beam, Particle* ptcle)
{
	int i = blockDim.x * blockIdx.x + threadIdx.x;
	
	const double beta0 = beam[_beta];

	const double ds    = drift[_length];
	const double ax    = drift[_apertureX];
	const double ay    = drift[_apertureY];
	
    if(ptcle[_f(i)])
	{
		double x0  = ptcle[_x(i)];
		double px0 = ptcle[_px(i)];
		double y0  = ptcle[_y(i)];
		double py0 = ptcle[_py(i)];
		double dp0 = ptcle[_dp(i)];
		
		double d1  = sqrt(1 - px0*px0 - py0*py0 + 2*dp0/beta0 + dp0*dp0);

		double x1  = x0 + ds*px0/d1;
		double y1  = y0 + ds*py0/d1;

		ptcle[_x(i)]   = x1;
		ptcle[_y(i)]   = y1;
		ptcle[_ct(i)] += ds*(1 - (1 + beta0*dp0)/d1)/beta0;
	    
		ptcle[_s(i)]  += ds;
		
		if( (ax>0) & (ay>0) )
			if((x1/ax)*(x1/ax) + (y1/ay)*(y1/ay) > 1)
				ptcle[_f(i)] = 0;
    }

	if(i==0)
		beam[_globaltime] += ds / (beta0*SpeedOfLight);
}


extern "C" __host__ void TrackDriftCUDA(DriftParameters_t drift)
{
	size_t parametersListSize = _nparams*sizeof(Component);

	Component* d_ComponentParameters;

	Component parameters[_nparams];
	parameters[_length]    = drift.length;
	parameters[_apertureX] = drift.apertureX;
	parameters[_apertureY] = drift.apertureY;

	cudaMalloc((void**)&d_ComponentParameters, parametersListSize);
	cudaMemcpy(d_ComponentParameters, parameters, parametersListSize, cudaMemcpyHostToDevice);

	DriftMap<<<blocksPerGrid, threadsPerBlock>>>(d_ComponentParameters, d_BeamParameters, d_ParticleArray);

	cudaFree(d_ComponentParameters);
}


extern "C" __host__ void CopyDriftCUDA(DriftParameters_t drift)
{
	size_t parametersListSize = _nparams*sizeof(Component);

	Component* d_ComponentParameters;

	Component parameters[_nparams];
	parameters[_length]    = drift.length;
	parameters[_apertureX] = drift.apertureX;
	parameters[_apertureY] = drift.apertureY;

	cudaMalloc((void**)&d_ComponentParameters, parametersListSize);
	cudaMemcpy(d_ComponentParameters, parameters, parametersListSize, cudaMemcpyHostToDevice);

	AppendComponent(&DriftMap, d_ComponentParameters);
}
