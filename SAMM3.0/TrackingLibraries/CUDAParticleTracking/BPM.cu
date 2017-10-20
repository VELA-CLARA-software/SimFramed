#include <cuda_runtime.h>

#include "CUDAParticleTracking.h"
#include "TrackCUDA.h"

extern "C" BeamParameters_t* beamParams;
extern Particle* d_ParticleArray;

double* d_bpmX;
double* d_bpmY;
int* d_totf;

inline int Div2(int n)
{
	return (n/2)+(n%2);
}


__global__ void BPMMap1(double* bpmX, double* bpmY, int* totf, Particle* ptcle, int imax, double ax, double ay)
{
	int i = blockDim.x * blockIdx.x + threadIdx.x;
	
	if(2*i<imax-1)
	{
		bpmX[i] = (ptcle[_x(2*i)] * ptcle[_f(2*i)]) + (ptcle[_x(2*i+1)] * ptcle[_f(2*i+1)]);
		bpmY[i] = (ptcle[_y(2*i)] * ptcle[_f(2*i)]) + (ptcle[_y(2*i+1)] * ptcle[_f(2*i+1)]);
		totf[i] = ptcle[_f(2*i)] + ptcle[_f(2*i+1)];
		
		if( (ax>0) & (ay>0) )
		{
			double x1 = ptcle[_x(2*i)];
			double y1 = ptcle[_y(2*i)];		
			if((x1/ax)*(x1/ax) + (y1/ay)*(y1/ay) > 1)
				ptcle[_f(2*i)] = 0;
				
			x1 = ptcle[_x(2*i+1)];
			y1 = ptcle[_y(2*i+1)];		
			if((x1/ax)*(x1/ax) + (y1/ay)*(y1/ay) > 1)
				ptcle[_f(2*i+1)] = 0;
		}
	}
		
	if(2*i==imax-1)
	{
		bpmX[i] = ptcle[_x(2*i)] * ptcle[_f(2*i)];
		bpmY[i] = ptcle[_y(2*i)] * ptcle[_f(2*i)];
		totf[i] = ptcle[_f(2*i)];
		
		if( (ax>0) & (ay>0) )
		{
			double x1 = ptcle[_x(2*i)];
			double y1 = ptcle[_y(2*i)];		
			if((x1/ax)*(x1/ax) + (y1/ay)*(y1/ay) > 1)
				ptcle[_f(2*i)] = 0;
		}
	}
}


__global__ void BPMMapN(double* bpmX, double* bpmY, int* totf, int n, int imax)
{
	int i = blockDim.x * blockIdx.x + threadIdx.x;
	
	if(2*n*i<imax-n)
	{
		bpmX[2*n*i] = bpmX[2*n*i] + bpmX[2*n*i + n];
		bpmY[2*n*i] = bpmY[2*n*i] + bpmY[2*n*i + n];
		totf[2*n*i] = totf[2*n*i] + totf[2*n*i + n];
	}
}


extern "C" __host__ void TrackBPMCUDA(BPMParameters_t bpm, double* bpmX, double* bpmY, int* nparticles)
{
	int imax    = beamParams->nparticles;
	int ntotals = Div2(imax);

	size_t sizeptcledouble = ntotals * sizeof(double);
	size_t sizeptcleint = ntotals * sizeof(int);

	cudaMalloc((void**)&d_bpmX, sizeptcledouble);
	cudaMalloc((void**)&d_bpmY, sizeptcledouble);
	cudaMalloc((void**)&d_totf, sizeptcleint);

	int blocksPerGrid = (ntotals + threadsPerBlock - 1) / threadsPerBlock;
	
	BPMMap1<<<blocksPerGrid, threadsPerBlock>>>(d_bpmX, d_bpmY, d_totf, d_ParticleArray, imax, bpm.apertureX, bpm.apertureY);
	
	imax = ntotals;
	int n = 1;
	while(ntotals>1)
	{
		ntotals = Div2(ntotals);
		blocksPerGrid = (ntotals + threadsPerBlock - 1) / threadsPerBlock;
		
		BPMMapN<<<blocksPerGrid, threadsPerBlock>>>(d_bpmX, d_bpmY, d_totf, n, imax);
		
		n = n * 2;
	}
	
	cudaMemcpy(bpmX, d_bpmX, sizeof(double), cudaMemcpyDeviceToHost);
	cudaMemcpy(bpmY, d_bpmY, sizeof(double), cudaMemcpyDeviceToHost);
	cudaMemcpy(nparticles, d_totf, sizeof(int), cudaMemcpyDeviceToHost);
	
	cudaFree(d_bpmX);
	cudaFree(d_bpmY);
	cudaFree(d_totf);
	
	if(*nparticles>0)
	{
		*bpmX /= *nparticles;
		*bpmY /= *nparticles;
	}

}
