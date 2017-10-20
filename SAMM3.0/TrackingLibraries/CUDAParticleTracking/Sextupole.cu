#include <cuda_runtime.h>
#include <math.h>

#include "CUDAParticleTracking.h"
#include "TrackCUDA.h"

#define _nparams 4
#define _length 0
#define _gradient 1
#define _apertureX 2
#define _apertureY 3

extern Beam* d_BeamParameters;
extern Particle* d_ParticleArray;
extern int blocksPerGrid;

__global__ void SextupoleMap(Component* sextupole, Beam* beam, Particle* ptcle)
{
	int i = blockDim.x * blockIdx.x + threadIdx.x;
	
	const double brho  = beam[_rigidity];
	const double beta0 = beam[_beta];

	const double ds = sextupole[_length];
	const double k2 = sextupole[_gradient] / brho;
	const double ax = sextupole[_apertureX];
	const double ay = sextupole[_apertureY];
	
    if(ptcle[_f(i)])
	{
		double px0 = ptcle[_px(i)];
		double py0 = ptcle[_py(i)];
		double dp0 = ptcle[_dp(i)];
		
		// First, apply a drift through ds/2
		double d1  = sqrt(1 - px0*px0 - py0*py0 + 2*dp0/beta0 + dp0*dp0);
		
		ptcle[_x(i)]  += ds*px0/d1/2;
		ptcle[_y(i)]  += ds*py0/d1/2;
		ptcle[_ct(i)] += ds*(1 - (1 + beta0*dp0)/d1)/beta0/2;
	    
	    // Then, apply a sextupole kick
	    double x1 = ptcle[_x(i)];
	    double y1 = ptcle[_y(i)];
        px0 += -(x1*x1 - y1*y1)*k2*ds/2;
        py0 += x1*y1*k2*ds;
	    
	    // Next, apply another drift through ds/2
		d1  = sqrt(1 - px0*px0 - py0*py0 + 2*dp0/beta0 + dp0*dp0);
		
		ptcle[_x(i)]  += ds*px0/d1/2;
		ptcle[_px(i)]  = px0;
		ptcle[_y(i)]  += ds*py0/d1/2;
		ptcle[_py(i)]  = py0;
		ptcle[_ct(i)] += ds*(1 - (1 + beta0*dp0)/d1)/beta0/2;
	    
	    // Finally, collimate
		ptcle[_s(i)] += ds;
		
		x1 = ptcle[_x(i)];
		y1 = ptcle[_y(i)];
		
		if( (ax>0) & (ay>0) )
			if((x1/ax)*(x1/ax) + (y1/ay)*(y1/ay) > 1)
				ptcle[_f(i)] = 0;
    }

	if(i==0)
		beam[_globaltime] += ds / (beta0 * SpeedOfLight);
}

extern "C" __host__ void TrackSextupoleCUDA(SextupoleParameters_t sextupole)
{
	size_t cptParamsListSize = _nparams*sizeof(Component);

	Component* d_ComponentParameters;

	Component parameters[_nparams];
	parameters[_length]    = sextupole.length;
	parameters[_gradient]  = sextupole.gradient;
	parameters[_apertureX] = sextupole.apertureX;
	parameters[_apertureY] = sextupole.apertureY;
		
	cudaMalloc((void**)&d_ComponentParameters, cptParamsListSize);

	cudaMemcpy(d_ComponentParameters, parameters, cptParamsListSize, cudaMemcpyHostToDevice);

	SextupoleMap<<<blocksPerGrid, threadsPerBlock>>>(d_ComponentParameters, d_BeamParameters, d_ParticleArray);

	cudaFree(d_ComponentParameters);
}
