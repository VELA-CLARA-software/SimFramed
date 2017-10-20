#include <cuda_runtime.h>
#include <math.h>

#include "CUDAParticleTracking.h"
#include "TrackCUDA.h"

#define _nparams 5
#define _length 0
#define _fieldX 1
#define _fieldY 2
#define _apertureX 3
#define _apertureY 4

extern Beam* d_BeamParameters;
extern Particle* d_ParticleArray;
extern int blocksPerGrid;

__global__ void OrbitCorrectorMap(Component* orbitcorrector, Beam* beam, Particle* ptcle)
{
	int i = blockDim.x * blockIdx.x + threadIdx.x;
	
	const double brho  = beam[_rigidity];
	const double beta0 = beam[_beta];

	const double ds    = orbitcorrector[_length];
	const double kx    = orbitcorrector[_fieldX] / brho;
	const double ky    = orbitcorrector[_fieldY] / brho;
	const double ax    = orbitcorrector[_apertureX];
	const double ay    = orbitcorrector[_apertureY];
	
	
    if(ptcle[_f(i)])
	{
		double x0  = ptcle[_x(i)];
		double px0 = ptcle[_px(i)];
		double y0  = ptcle[_y(i)];
		double py0 = ptcle[_py(i)];
		double dp0 = ptcle[_dp(i)];
		double ct0 = ptcle[_ct(i)];
		
		double d1  = sqrt(1 + 2*dp0/beta0 + dp0*dp0);

		double x1  = x0 + ds*px0/d1 - ds*ds*ky/d1/2;
		double px1 =         px0    - ds*ky;

		double y1  = y0 + ds*py0/d1 + ds*ds*kx/d1/2;
		double py1 =         py0    + ds*kx;

		double f1  = ds*(1/beta0 + dp0)/d1/d1/d1/2;

		double c0  = ct0 + ds/beta0 - ds*(1/beta0 + dp0)/d1 -
				      ds*ds*f1*(kx*kx + ky*ky)/3;

		double ct1 = c0 + ds*f1*(ky*px0 - kx*py0) - f1*(px0*px0 + py0*py0);

		ptcle[_x(i)]  = x1;
		ptcle[_px(i)] = px1;
		ptcle[_y(i)]  = y1;
		ptcle[_py(i)] = py1;
		ptcle[_ct(i)] = ct1;
	    
	    // Finally, collimate
		ptcle[_s(i)] += ds;
		
		if( (ax>0) & (ay>0) )
			if((x1/ax)*(x1/ax) + (y1/ay)*(y1/ay) > 1)
				ptcle[_f(i)] = 0;
    }

	if(i==0)
		beam[_globaltime] += ds / (beta0 * SpeedOfLight);
}


extern "C" __host__ void TrackOrbitCorrectorCUDA(OrbitCorrectorParameters_t orbitcorrector)
{
	size_t cptParamsListSize = _nparams*sizeof(Component);

	Component* d_ComponentParameters;

	Component parameters[_nparams];
	parameters[_length]    = orbitcorrector.length;
	parameters[_fieldX]    = orbitcorrector.fieldX;
	parameters[_fieldY]    = orbitcorrector.fieldY;
	parameters[_apertureX] = orbitcorrector.apertureX;
	parameters[_apertureY] = orbitcorrector.apertureY;

	cudaMalloc((void**)&d_ComponentParameters, cptParamsListSize);
	cudaMemcpy(d_ComponentParameters, parameters, cptParamsListSize, cudaMemcpyHostToDevice);

	OrbitCorrectorMap<<<blocksPerGrid, threadsPerBlock>>>(d_ComponentParameters, d_BeamParameters, d_ParticleArray);

	cudaFree(d_ComponentParameters);
}


extern "C" __host__ void CopyOrbitCorrectorCUDA(OrbitCorrectorParameters_t orbitcorrector)
{
	size_t cptParamsListSize = _nparams*sizeof(Component);

	Component* d_ComponentParameters;

	Component parameters[_nparams];
	parameters[_length]    = orbitcorrector.length;
	parameters[_fieldX]    = orbitcorrector.fieldX;
	parameters[_fieldY]    = orbitcorrector.fieldY;
	parameters[_apertureX] = orbitcorrector.apertureX;
	parameters[_apertureY] = orbitcorrector.apertureY;

	cudaMalloc((void**)&d_ComponentParameters, cptParamsListSize);
	cudaMemcpy(d_ComponentParameters, parameters, cptParamsListSize, cudaMemcpyHostToDevice);

	AppendComponent(&OrbitCorrectorMap, d_ComponentParameters);
}