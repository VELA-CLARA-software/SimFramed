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

__global__ void QuadrupoleMap(Component* quadrupole, Beam* beam, Particle* ptcle)
{
	int i = blockDim.x * blockIdx.x + threadIdx.x;
	
	double beta0 = beam[_beta];
	double brho  = beam[_rigidity];

	double ds = quadrupole[_length];
	double k1 = quadrupole[_gradient] / brho;
	
	double ax = quadrupole[_apertureX];
	double ay = quadrupole[_apertureY];
	
    if(ptcle[_f(i)])
	{
		double x0  = ptcle[_x(i)];
		double px0 = ptcle[_px(i)];
		double y0  = ptcle[_y(i)];
		double py0 = ptcle[_py(i)];
		double ct0 = ptcle[_ct(i)];
		double dp0 = ptcle[_dp(i)];
	
		double d1  = sqrt(1 + 2*dp0/beta0 + dp0*dp0);

		double w   = sqrt(fabs(k1)/d1);

		double xs;
		double xc;
		double ys;
		double yc;
		double xs2;
		double ys2;
	    
		if(k1>=0)
		{
		   xs  = sin(w*ds);
		   xc  = cos(w*ds);
		   ys  = sinh(w*ds);
		   yc  = cosh(w*ds);
		   xs2 = sin(2*w*ds);
		   ys2 = sinh(2*w*ds);
		}
		else
		{
		   xs  = sinh(w*ds);
		   xc  = cosh(w*ds);
		   ys  = sin(w*ds);
		   yc  = cos(w*ds);
		   xs2 = sinh(2*w*ds);
		   ys2 = sin(2*w*ds);
		}

		double x1  =  x0*xc      + px0*xs*w/fabs(k1);
		double px1 = -k1*x0*xs/w + px0*xc;
		double y1  =  y0*yc      + py0*ys*w/fabs(k1);
		double py1 =  k1*y0*ys/w + py0*yc;
		
		double d0  =  1/beta0 + dp0;
		double d2  = -d0/d1/d1/d1/2;
	    
	    double c0  =  (1/beta0 - d0/d1)*ds;
		double c11 =  k1*k1*d2*(xs2/w - 2*ds)/w/w/4;
		double c12 = -k1*d2*xs*xs/w/w;
		double c22 =  d2*(xs2/w + 2*ds)/4;
		double c33 =  k1*k1*d2*(ys2/w - 2*ds)/w/w/4;
		double c34 =  k1*d2*ys*ys/w/w;
		double c44 =  d2*(ys2/w + 2*ds)/4;
	    
		double ct1 = ct0 + c0;
		double ctx = c11*x0*x0 + c12*x0*px0 + c22*px0*px0;
		double cty = c33*y0*y0 + c34*y0*py0 + c44*py0*py0;
						 
		ptcle[_x(i)]  =  x1;
		ptcle[_px(i)] = px1;	
		ptcle[_y(i)]  =  y1;
		ptcle[_py(i)] = py1;
		ptcle[_ct(i)] = ct1 + ctx + cty;
	    
		ptcle[_s(i)] += ds;
		
		if( (ax>0) & (ay>0) )
			if((x1/ax)*(x1/ax) + (y1/ay)*(y1/ay) > 1)
				ptcle[_f(i)] = 0;
	}

	if(i==0)
		beam[_globaltime] += ds / (beta0 * SpeedOfLight);
}

extern "C" __host__ void TrackQuadrupoleCUDA(QuadrupoleParameters_t quadrupole)
{
	size_t cptParamsListSize = _nparams*sizeof(Component);

	Component* d_ComponentParameters;

	Component parameters[_nparams];
	parameters[_length]    = quadrupole.length;
	parameters[_gradient]  = quadrupole.gradient;
	parameters[_apertureX] = quadrupole.apertureX;
	parameters[_apertureY] = quadrupole.apertureY;

	cudaMalloc((void**)&d_ComponentParameters, cptParamsListSize);

	cudaMemcpy(d_ComponentParameters, parameters, cptParamsListSize, cudaMemcpyHostToDevice);

	QuadrupoleMap<<<blocksPerGrid, threadsPerBlock>>>(d_ComponentParameters, d_BeamParameters, d_ParticleArray);

	cudaFree(d_ComponentParameters);	
}