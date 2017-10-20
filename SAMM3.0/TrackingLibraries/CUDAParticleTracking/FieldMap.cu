#include <cuda_runtime.h>
#include <math.h>
#include <stdlib.h>

#include "CUDAParticleTracking.h"
#include "TrackCUDA.h"

extern "C" BeamParameters_t* beamParams;

extern Beam* d_BeamParameters;
extern Particle* d_ParticleArray;
extern int blocksPerGrid;

#define interp3(B) \
	(a000*B[(nx-1) + sx*(ny-1) + sx*sy*(nz-1)] + \
	 a001*B[(nx-1) + sx*(ny-1) + sx*sy* nz   ] + \
	 a010*B[(nx-1) + sx* ny    + sx*sy*(nz-1)] + \
	 a011*B[(nx-1) + sx* ny    + sx*sy* nz   ] + \
	 a100*B[ nx    + sx*(ny-1) + sx*sy*(nz-1)] + \
	 a101*B[ nx    + sx*(ny-1) + sx*sy* nz   ] + \
	 a110*B[ nx    + sx* ny    + sx*sy*(nz-1)] + \
	 a111*B[ nx    + sx* ny    + sx*sy* nz   ])

__global__ void FieldMaptcletep(Component* fieldmap, double* gridX, double* gridY, double* gridZ, double* Bx, double* By, double* Bz, Particle* ptcle)
{
	int i = blockDim.x * blockIdx.x + threadIdx.x;
	
	const double length = fieldmap[0];
	const double ax     = fieldmap[1];
	const double ay     = fieldmap[2];
	
	const double gamma0 = fieldmap[3];
	const double qdmc   = fieldmap[4];
	
	const int nsteptcle   = (int)fieldmap[5];
	
	const int sx       = (int)fieldmap[6];
	const int sy       = (int)fieldmap[7];
	const int sz       = (int)fieldmap[8];
	
	const double b0g0   = sqrt(gamma0*gamma0 - 1);
	const double beta0  = b0g0/gamma0;
	const double bmin   = 1e-9;
	
    if(ptcle[_f(i)])
	{
		double x0    = ptcle[_x(i)];
		double px0   = ptcle[_px(i)];
		double y0    = ptcle[_y(i)];
		double py0   = ptcle[_py(i)];
		double ct0   = ptcle[_ct(i)];
		double dp0   = ptcle[_dp(i)];
		
		double z0    = 0;
		
		double gamma = gamma0 + b0g0*dp0;
		double bx0   = b0g0*px0/gamma;
		double by0   = b0g0*py0/gamma;
		double bz0   = sqrt(1 - bx0*bx0 - by0*by0 - 1/gamma/gamma);
		
		double k     = qdmc/gamma;
		
		int   nx    = (int)(sx/2);
		int   ny    = (int)(sy/2);
		int   nz    = (int)(sz/2);
		
		double s1    = 0;
		double s2    = 0;
		double c1    = 0;
		double c2    = 0;
		
		for(int n=0; n<nsteptcle; n++)
		{
			double cdt  = (length - z0)/(nsteptcle - n)/bz0;

			if(x0>gridX[nx])
				while(x0>gridX[nx]   & nx<sx) nx++;
			else
				while(x0<gridX[nx-1] & nx>1)  nx--;
				
			if(y0>gridY[ny])
				while(y0>gridY[ny]   & ny<sy) ny++;
			else
				while(y0<gridY[ny-1] & ny>1)  ny--;
			
			if(z0>gridZ[nz])
				while(z0>gridZ[nz]   & nz<sz) nz++;
			else
				while(z0<gridZ[nz-1] & nz>1)  nz--;

			__syncthreads();
			
			double ax1  = (x0 - gridX[nx-1])/(gridX[nx] - gridX[nx-1]);
			double ay1  = (y0 - gridY[ny-1])/(gridY[ny] - gridY[ny-1]);
			double az1  = (z0 - gridZ[nz-1])/(gridZ[nz] - gridZ[nz-1]);
			
			double a000 = (1-ax1)*(1-ay1)*(1-az1);
			double a001 = (1-ax1)*(1-ay1)*   az1 ;
			double a010 = (1-ax1)*   ay1 *(1-az1);
			double a011 = (1-ax1)*   ay1 *   az1 ;
			double a100 =    ax1 *(1-ay1)*(1-az1);
			double a101 =    ax1 *(1-ay1)*   az1 ;
			double a110 =    ax1 *   ay1 *(1-az1);
			double a111 =    ax1 *   ay1 *   az1 ;
			
			double Bx0  = interp3(Bx)*k;
			double By0  = interp3(By)*k;
			double Bz0  = interp3(Bz)*k;
			
			double bdotv = Bx0*bx0 + By0*by0 + Bz0*bz0;
			double bmag  = sqrt(Bx0*Bx0 + By0*By0 + Bz0*Bz0);
			
			if(abs(bmag)<bmin)
			{
				s1		= cdt;
				s2      = bdotv*cdt*cdt*cdt/6;
				c1      = 1;
				c2      = cdt*cdt/2;			
			}
			else
			{			
				s1      = sin(bmag*cdt)/bmag;
				s2      = bdotv*(bmag*cdt - sin(bmag*cdt))/bmag/bmag/bmag;
				c1      = cos(bmag*cdt);
				c2      = (1 - c1)/bmag/bmag;
			}
			
            double x1    = x0 + bx0*s1 + (by0*Bz0 - bz0*By0)*c2 + Bx0*s2;
            double y1    = y0 + by0*s1 + (bz0*Bx0 - bx0*Bz0)*c2 + By0*s2;
            double z1    = z0 + bz0*s1 + (bx0*By0 - by0*Bx0)*c2 + Bz0*s2;
            
            double bx1   = bx0*c1 + (by0*Bz0 - bz0*By0)*s1 + Bx0*bdotv*c2;
            double by1   = by0*c1 + (bz0*Bx0 - bx0*Bz0)*s1 + By0*bdotv*c2;
            double bz1   = bz0*c1 + (bx0*By0 - by0*Bx0)*s1 + Bz0*bdotv*c2;

            ct0         = ct0 + (bz0/beta0 - 1)*cdt;

            x0  = x1;
            y0  = y1;
            z0  = z1;

            bx0 = bx1;
            by0 = by1;
            bz0 = bz1;
		}
		
		ptcle[_x(i)]  = x0;
		ptcle[_px(i)] = bx0*gamma/b0g0;
		ptcle[_y(i)]  = y0;
		ptcle[_py(i)] = by0*gamma/b0g0;
		ptcle[_ct(i)] = ct0;
	    
		ptcle[_s(i)] += length;
		
		if( (ax>0) & (ay>0) )
			if((x0/ax)*(x0/ax) + (y0/ay)*(y0/ay) > 1)
				ptcle[_f(i)] = 0;
    }
}


#define CopyGrid(gsize, grid, d_grid) \
	sizeGrid = gsize * sizeof(double); \
	h_grid = (double*)malloc(sizeGrid); \
	for(int n=0; n<gsize; n++) h_grid[n] = (double)grid[n]; \
	cudaMalloc((void**)&d_grid, sizeGrid); \
	cudaMemcpy(d_grid, h_grid, sizeGrid, cudaMemcpyHostToDevice); \
	free(h_grid);


extern "C" __host__ void TrackFieldMapCUDA(FieldMapParameters_t fieldmap)
{
	size_t sizeGrid;
	double* h_grid;

	double* d_gridX;
	double* d_gridY;
	double* d_gridZ;

	CopyGrid(fieldmap.gridXsize, fieldmap.gridX, d_gridX);
	CopyGrid(fieldmap.gridYsize, fieldmap.gridY, d_gridY);
	CopyGrid(fieldmap.gridZsize, fieldmap.gridZ, d_gridZ);
	
	int fieldsize = fieldmap.gridXsize * fieldmap.gridYsize * fieldmap.gridZsize;

	double* d_Bx;
	double* d_By;
	double* d_Bz;

	CopyGrid(fieldsize, fieldmap.Bx, d_Bx);
	CopyGrid(fieldsize, fieldmap.By, d_By);
	CopyGrid(fieldsize, fieldmap.Bz, d_Bz);
	
	double brho   = beamParams->rigidity;
	double qdmc   = beamParams->charge / beamParams->mass / SpeedOfLight;
	double bg     = brho * qdmc;
	double gamma0 = sqrt(1+bg*bg);

	size_t sizeParams = 9*sizeof(Component);
	Component* fieldmapParams = (Component*)malloc(sizeParams);
	fieldmapParams[0] = (Component)fieldmap.length;
	fieldmapParams[1] = (Component)fieldmap.apertureX;
	fieldmapParams[2] = (Component)fieldmap.apertureY;
	fieldmapParams[3] = (Component)gamma0;
	fieldmapParams[4] = (Component)qdmc;
	fieldmapParams[5] = (Component)fieldmap.nsteptcle;
	fieldmapParams[6] = (Component)fieldmap.gridXsize;
	fieldmapParams[7] = (Component)fieldmap.gridYsize;
	fieldmapParams[8] = (Component)fieldmap.gridZsize;
		
	Component* d_Params;
	cudaMalloc((void**)&d_Params, sizeParams);
	cudaMemcpy(d_Params, fieldmapParams, sizeParams, cudaMemcpyHostToDevice);

	// Invoke kernel
	FieldMaptcletep<<<blocksPerGrid, threadsPerBlock>>>(d_Params, d_gridX, d_gridY, d_gridZ, d_Bx, d_By, d_Bz, d_ParticleArray);
	
	double beta0  = sqrt(1 - 1/gamma0/gamma0);
	beamParams->globaltime += fieldmap.length / (beta0*SpeedOfLight);

	free(fieldmapParams);
	
	cudaFree(d_Params);
	
	cudaFree(d_gridX);
	cudaFree(d_gridY);
	cudaFree(d_gridZ);
	
	cudaFree(d_Bx);
	cudaFree(d_By);
	cudaFree(d_Bz);
}
