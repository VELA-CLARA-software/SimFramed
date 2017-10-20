#include <cuda_runtime.h>
#include <math.h>

#include "CUDAParticleTracking.h"
#include "TrackCUDA.h"

#define _nparams 5
#define _length 0
#define _field 1
#define _taper 2
#define _apertureX 3
#define _apertureY 4

extern Beam* d_BeamParameters;
extern Particle* d_ParticleArray;
extern int blocksPerGrid;

__global__ void SolenoidMap(Component* solenoid, Beam* beam, Particle* ptcle)
{
	int i = blockDim.x * blockIdx.x + threadIdx.x;
	
	const double ds    = solenoid[_length];
	const double b0    = solenoid[_field];
	const double g     = solenoid[_taper];
	
	const double ax    = solenoid[_apertureX];
	const double ay    = solenoid[_apertureY];
	
	const double brho  = beam[_rigidity];
	const double beta0 = beam[_beta];

    if(ptcle[_f(i)])
	{
		double x0  = ptcle[_x(i)];
		double px0 = ptcle[_px(i)];
		double y0  = ptcle[_y(i)];
		double py0 = ptcle[_py(i)];
		double dp0 = ptcle[_dp(i)];
		double ct0 = ptcle[_ct(i)];
		
		
        const double gds = g*ds;
        const double b1  = b0 / (1 + gds);
            
        const double helmA = sqrt(b0/b1);
        const double helmB = 2*brho*(1+dp0)/sqrt(b0*b1);
        const double helmF = -1 / helmB;
        const double helmG =  1 / helmA;
            
        double w = b0*ds/2/brho/(1+dp0);
        if(gds!=0)
            w = w*log(1+gds)/gds;
            
        double cw2 = cos(w)*cos(w);
        double s2w = sin(2*w);
        double sw2 = sin(w)*sin(w);
        
		double d1  = sqrt(1 - px0*px0 - py0*py0
		                  + 2*dp0/beta0 + dp0*dp0);
		                  
		double ct1 = ct0 + ds*(1 - (1+beta0*dp0)/d1)/beta0/2;
		                  
        double x1  = helmA*cw2*x0   + helmB*s2w*px0/2 + helmA*s2w*y0/2 + helmB*sw2*py0;
        double px1 = helmF*s2w*x0/2 + helmG*cw2*px0   + helmF*sw2*y0   + helmG*s2w*py0/2;
        double y1  =-helmA*s2w*x0/2 - helmB*sw2*px0   + helmA*cw2*y0   + helmB*s2w*py0/2;
        double py1 =-helmF*sw2*x0   - helmG*s2w*px0/2 + helmF*s2w*y0/2 + helmG*cw2*py0;

		d1  = sqrt(1 - px0*px0 - py0*py0
		             + 2*dp0/beta0 + dp0*dp0);
		                  
		ct1 += ds*(1 - (1+beta0*dp0)/d1)/beta0/2;

		ptcle[_x(i)]  = x1;
		ptcle[_px(i)] = px1;	
		ptcle[_y(i)]  = y1;
		ptcle[_py(i)] = py1;
		ptcle[_ct(i)] = ct1;
	    
		ptcle[_s(i)] += ds;
		
		if( (ax>0) & (ay>0) )
			if((x1/ax)*(x1/ax) + (y1/ay)*(y1/ay) > 1)
				ptcle[_f(i)] = 0;
	}

	if(i==0)
		beam[_globaltime] += ds / (beta0 * SpeedOfLight);
}


extern "C" __host__ void TrackSolenoidCUDA(SolenoidParameters_t solenoid)
{
	size_t cptParamsListSize = _nparams*sizeof(Component);

	Component* d_ComponentParameters;

	Component parameters[_nparams];
	parameters[_length]    = solenoid.length;
	parameters[_field]     = solenoid.field;
	parameters[_taper]     = solenoid.taper;
	parameters[_apertureX] = solenoid.apertureX;
	parameters[_apertureY] = solenoid.apertureY;

	cudaMalloc((void**)&d_ComponentParameters, cptParamsListSize);

	cudaMemcpy(d_ComponentParameters, parameters, cptParamsListSize, cudaMemcpyHostToDevice);

	SolenoidMap<<<blocksPerGrid, threadsPerBlock>>>(d_ComponentParameters, d_BeamParameters, d_ParticleArray);

	cudaFree(d_ComponentParameters);
}
