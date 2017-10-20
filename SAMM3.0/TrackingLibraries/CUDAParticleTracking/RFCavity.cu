#include <cuda_runtime.h>
#include <math.h>

#include "CUDAParticleTracking.h"
#include "TrackCUDA.h"

#define pi (double)3.1415926536

#define _nparams 7
#define _length 0
#define _voltage 1
#define _frequency 2
#define _phase 3
#define _apertureX 4
#define _apertureY 5
#define _masteroscillatorfrequency 6

extern Beam* d_BeamParameters;
extern Particle* d_ParticleArray;
extern int blocksPerGrid;

__global__ void RFCavityMap(Component* rfcavity, Beam* beam, Particle* ptcle)
{
	int i = blockDim.x * blockIdx.x + threadIdx.x;

	const double brho   = beam[_rigidity];
	const double beta0  = beam[_beta];
	double globaltime   = beam[_globaltime];
	
	const double ds     = rfcavity[_length];
	const double v0     = rfcavity[_voltage] / brho / SpeedOfLight;
	const double f      = rfcavity[_frequency];
	const double phi    = rfcavity[_phase];

	const double ax     = rfcavity[_apertureX];
	const double ay     = rfcavity[_apertureY];
	

    if(ptcle[_f(i)])
	{
		const double mofreq = rfcavity[_masteroscillatorfrequency];
		int p = (int)(globaltime * mofreq);
		globaltime -= (double)(p) / mofreq;

		double x0  = ptcle[_x(i)];
		double px0 = ptcle[_px(i)];
		double y0  = ptcle[_y(i)];
		double py0 = ptcle[_py(i)];
		double ct0 = ptcle[_ct(i)];
		double dp0 = ptcle[_dp(i)];

		// First, apply a drift through ds/2
		double d1  = sqrt(1 - px0*px0 - py0*py0 + 2*dp0/beta0 + dp0*dp0);

		double x1  = x0  + ds*px0/d1/2;
		double y1  = y0  + ds*py0/d1/2;
		double ct1 = ct0 + ds*(1 - (1 + beta0*dp0)/d1)/beta0/2;
	    
		// Then, apply an RF kick
		double t  = globaltime - ct1/beta0/SpeedOfLight;
        dp0      += v0*sin(2*pi*f*t + phi);

		// Next, apply another drift through ds/2
		d1  = sqrt(1 - px0*px0 - py0*py0 + 2*dp0/beta0 + dp0*dp0);

		x1  += ds*px0/d1/2;
		y1  += ds*py0/d1/2;
		ct1 += ds*(1 - (1 + beta0*dp0)/d1)/beta0/2;
		
		ptcle[_x(i)]  = x1;
		ptcle[_y(i)]  = y1;
		ptcle[_ct(i)] = ct1;
		ptcle[_dp(i)] = dp0;
		
		ptcle[_s(i)] += ds;
		
		// Finally, collimate
		if( (ax>0) & (ay>0) )
			if((x1/ax)*(x1/ax) + (y1/ay)*(y1/ay) > 1)
				ptcle[_f(i)] = 0;
    }

	if(i==0)
		beam[_globaltime] = globaltime + ds/beta0/SpeedOfLight;

}

extern "C" __host__ void TrackRFCavityCUDA(RFCavityParameters_t rfcavity)
{
	size_t cptParamsListSize = _nparams*sizeof(Component);

	Component* d_ComponentParameters;

	Component parameters[_nparams];
	parameters[_length]    = rfcavity.length;
	parameters[_voltage]   = rfcavity.voltage;
	parameters[_frequency] = rfcavity.frequency;
	parameters[_phase]     = rfcavity.phase;
	parameters[_apertureX] = rfcavity.apertureX;
	parameters[_apertureY] = rfcavity.apertureY;
	parameters[_masteroscillatorfrequency] = rfcavity.masteroscillatorfrequency;

	cudaMalloc((void**)&d_ComponentParameters, cptParamsListSize);
	cudaMemcpy(d_ComponentParameters, parameters, cptParamsListSize, cudaMemcpyHostToDevice);

	RFCavityMap<<<blocksPerGrid, threadsPerBlock>>>(d_ComponentParameters, d_BeamParameters, d_ParticleArray);

	cudaFree(d_ComponentParameters);
}

extern "C" __host__ void CopyRFCavityCUDA(RFCavityParameters_t rfcavity)
{
	size_t cptParamsListSize = _nparams*sizeof(Component);

	Component* d_ComponentParameters;

	Component parameters[_nparams];
	parameters[_length]    = rfcavity.length;
	parameters[_voltage]   = rfcavity.voltage;
	parameters[_frequency] = rfcavity.frequency;
	parameters[_phase]     = rfcavity.phase;
	parameters[_apertureX] = rfcavity.apertureX;
	parameters[_apertureY] = rfcavity.apertureY;
	parameters[_masteroscillatorfrequency] = rfcavity.masteroscillatorfrequency;

	cudaMalloc((void**)&d_ComponentParameters, cptParamsListSize);
	cudaMemcpy(d_ComponentParameters, parameters, cptParamsListSize, cudaMemcpyHostToDevice);
	
	AppendComponent(&RFCavityMap, d_ComponentParameters);
}
