#include <cuda_runtime.h>
#include <math.h>
#include <stdlib.h>

#include "CUDAParticleTracking.h"
#include "TrackCUDA.h"

#define _nparams 5
#define _angle 0
#define _curvature 1
#define _ncpts 2
#define _apertureX 3
#define _apertureY 4

extern Beam* d_BeamParameters;
extern Particle* d_ParticleArray;
extern int blocksPerGrid;

__global__ void MultipoleMap(Component* multipole, Beam* beam, Particle* ptcle)
{
	int i = blockDim.x * blockIdx.x + threadIdx.x;
	
	const double brho      = beam[_rigidity];

	const double bendangle =      multipole[_angle];
	const double curvature =      multipole[_curvature];
	const int    ncpts     = (int)multipole[_ncpts];
	const double ax        =      multipole[_apertureX];
	const double ay        =      multipole[_apertureY];


    if(ptcle[_f(i)])
	{
		double x0  = ptcle[_x(i)];
		double px0 = ptcle[_px(i)];
		double y0  = ptcle[_y(i)];
		double py0 = ptcle[_py(i)];
		double ct0 = ptcle[_ct(i)];
		double dp0 = ptcle[_dp(i)];

		double dpx = 0;
		double dpy = 0;

		double ptclei  = atan2(y0,x0);
		double r    = sqrt(x0*x0 + y0*y0);

		int    indx = 0;

		for(indx=0; indx<ncpts; indx++)
		{
			int findex = (int)multipole[_nparams + indx          ];
			double bnL =      multipole[_nparams + indx + ncpts  ];
			double anL =      multipole[_nparams + indx + ncpts*2];
			double cs  = cos(findex*ptclei);
			double sn  = sin(findex*ptclei);
			double pr  = pow(r,findex) / brho;
			dpx       += pr*(anL*sn - bnL*cs);
			dpy       += pr*(bnL*sn + anL*cs);
		}

		ptcle[_px(i)]  = px0 + dpx + bendangle*(1 + dp0 - curvature*x0);
		ptcle[_py(i)]  = py0 + dpy;
		ptcle[_ct(i)]  = ct0 - bendangle*x0;
	    
	    // Finally, collimate
		if( (ax>0) & (ay>0) )
			if((x0/ax)*(x0/ax) + (y0/ay)*(y0/ay) > 1)
				ptcle[_f(i)] = 0;
    }

}


extern "C" __host__ int factorial(int n)
{
	return (n==1 || n==0) ? 1 : factorial(n - 1) * n;
}


extern "C" __host__ void TrackMultipoleCUDA(MultipoleParameters_t multipole)
{
	int ncpts = multipole.ncpts;

	size_t cptParamsListSize = (_nparams + 3*ncpts) * sizeof(Component);

	Component* d_ComponentParameters;

	Component* parameters  = (Component*)malloc(cptParamsListSize);
	parameters[_angle]     = multipole.angle;
	parameters[_curvature] = multipole.curvature;
	parameters[_ncpts]     = (Component)ncpts;
	parameters[_apertureX] = multipole.apertureX;
	parameters[_apertureY] = multipole.apertureY;

	for(int indx=0; indx<ncpts; indx++)
	{
		int findex = (int)multipole.fieldindex[indx];
		int findexfact = factorial(findex);
		parameters[_nparams + indx          ] = (Component)findex;
		parameters[_nparams + indx + ncpts  ] = multipole.bnL[indx] / findexfact;
		parameters[_nparams + indx + ncpts*2] = multipole.anL[indx] / findexfact;
	}

	cudaMalloc((void**)&d_ComponentParameters, cptParamsListSize);
	cudaMemcpy(d_ComponentParameters, parameters, cptParamsListSize, cudaMemcpyHostToDevice);

	MultipoleMap<<<blocksPerGrid, threadsPerBlock>>>(d_ComponentParameters, d_BeamParameters, d_ParticleArray);

	cudaFree(d_ComponentParameters);
}


extern "C" __host__ void CopyMultipoleCUDA(MultipoleParameters_t multipole)
{
	int ncpts = multipole.ncpts;

	size_t cptParamsListSize = (_nparams + 3*ncpts) * sizeof(Component);

	Component* d_ComponentParameters;

	Component* parameters = (Component*)malloc(cptParamsListSize);
	parameters[_angle]     = multipole.angle;
	parameters[_curvature] = multipole.curvature;
	parameters[_ncpts]     = (Component)ncpts;
	parameters[_apertureX] = multipole.apertureX;
	parameters[_apertureY] = multipole.apertureY;

	for(int indx=0; indx<ncpts; indx++)
	{
		int findex = (int)multipole.fieldindex[indx];
		int findexfact = factorial(findex);
		parameters[_nparams + indx          ] = (Component)findex;
		parameters[_nparams + indx + ncpts  ] = multipole.bnL[indx] / findexfact;
		parameters[_nparams + indx + ncpts*2] = multipole.anL[indx] / findexfact;
	}

	cudaMalloc((void**)&d_ComponentParameters, cptParamsListSize);
	cudaMemcpy(d_ComponentParameters, parameters, cptParamsListSize, cudaMemcpyHostToDevice);

	AppendComponent(&MultipoleMap, d_ComponentParameters);
}
