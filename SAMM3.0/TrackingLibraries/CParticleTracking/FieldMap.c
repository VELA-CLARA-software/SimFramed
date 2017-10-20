
#include <math.h>

#include "CParticleTracking.h"
#include "GlobalVariables.h"

#define interp3(B) \
	(a000*B[(nx-1) + sx*(ny-1) + sx*sy*(nz-1)] + \
	 a001*B[(nx-1) + sx*(ny-1) + sx*sy* nz   ] + \
	 a010*B[(nx-1) + sx* ny    + sx*sy*(nz-1)] + \
	 a011*B[(nx-1) + sx* ny    + sx*sy* nz   ] + \
	 a100*B[ nx    + sx*(ny-1) + sx*sy*(nz-1)] + \
	 a101*B[ nx    + sx*(ny-1) + sx*sy* nz   ] + \
	 a110*B[ nx    + sx* ny    + sx*sy*(nz-1)] + \
	 a111*B[ nx    + sx* ny    + sx*sy* nz   ])

void TrackFieldMap(FieldMapParameters_t fieldmap)
{
	const double ds = fieldmap.length;
	const double ax = fieldmap.apertureX;
	const double ay = fieldmap.apertureY;

	const double length = fieldmap.length;

	const double qdmc   = beamParams->charge / beamParams->mass / SpeedOfLight;
	const double b0g0   = beamParams->rigidity * qdmc;
	const double gamma0 = sqrt(1+b0g0*b0g0);
	const double beta0  = b0g0/gamma0;

	double* gridX = fieldmap.gridX;
	double* gridY = fieldmap.gridY;
	double* gridZ = fieldmap.gridZ;

	double* Bx    = fieldmap.Bx;
	double* By    = fieldmap.By;
	double* Bz    = fieldmap.Bz;

	const int sx  = fieldmap.gridXsize;
	const int sy  = fieldmap.gridYsize;
	const int sz  = fieldmap.gridZsize;

	const int nsteps = fieldmap.nsteps;

	const double bmin = 1e-9;

	int n;
	for(n=0; n<beamParams->nparticles; n++)
	{
		if(distance[n*2 + 1])
		{

			double x0  = beamPS[n*6];
			double px0 = beamPS[n*6 + 1];
			double y0  = beamPS[n*6 + 2];
			double py0 = beamPS[n*6 + 3];
			double ct0 = beamPS[n*6 + 4];
			double dp0 = beamPS[n*6 + 5];

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

			double ax1;
			double ay1;
			double az1;

			double a000;
			double a001;
			double a010;
			double a011;
			double a100;
			double a101;
			double a110;
			double a111;

			double Bx0;
			double By0;
			double Bz0;

			double bdotv;
			double bmag;

			double x1;
			double y1;
			double z1;

			double bx1;
			double by1;
			double bz1;

			double px1;
			double py1;

			double theta;
			double phi;

			int ns;

			for(ns=0; ns<nsteps; ns++)
			{
				double cdt  = (length - z0)/(nsteps - ns)/bz0;

				if(x0>gridX[nx])
					while(x0>gridX[nx]   && nx<sx) nx++;
				else
					while(x0<gridX[nx-1] && nx>1)  nx--;
					
				if(y0>gridY[ny])
					while(y0>gridY[ny]   && ny<sy) ny++;
				else
					while(y0<gridY[ny-1] && ny>1)  ny--;
				
				if(z0>gridZ[nz])
					while(z0>gridZ[nz]   && nz<sz) nz++;
				else
					while(z0<gridZ[nz-1] && nz>1)  nz--;

				ax1  = (x0 - gridX[nx-1])/(gridX[nx] - gridX[nx-1]);
				ay1  = (y0 - gridY[ny-1])/(gridY[ny] - gridY[ny-1]);
				az1  = (z0 - gridZ[nz-1])/(gridZ[nz] - gridZ[nz-1]);
				
				a000 = (1-ax1)*(1-ay1)*(1-az1);
				a001 = (1-ax1)*(1-ay1)*   az1 ;
				a010 = (1-ax1)*   ay1 *(1-az1);
				a011 = (1-ax1)*   ay1 *   az1 ;
				a100 =    ax1 *(1-ay1)*(1-az1);
				a101 =    ax1 *(1-ay1)*   az1 ;
				a110 =    ax1 *   ay1 *(1-az1);
				a111 =    ax1 *   ay1 *   az1 ;
				
				Bx0  = interp3(Bx)*k;
				By0  = interp3(By)*k;
				Bz0  = interp3(Bz)*k;
				
				bdotv = Bx0*bx0 + By0*by0 + Bz0*bz0;
				bmag  = sqrt(Bx0*Bx0 + By0*By0 + Bz0*Bz0);
				
				if(fabs(bmag)<bmin)
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
				
				x1    = x0 + bx0*s1 + (by0*Bz0 - bz0*By0)*c2 + Bx0*s2;
				y1    = y0 + by0*s1 + (bz0*Bx0 - bx0*Bz0)*c2 + By0*s2;
				z1    = z0 + bz0*s1 + (bx0*By0 - by0*Bx0)*c2 + Bz0*s2;
	            
				bx1   = bx0*c1 + (by0*Bz0 - bz0*By0)*s1 + Bx0*bdotv*c2;
				by1   = by0*c1 + (bz0*Bx0 - bx0*Bz0)*s1 + By0*bdotv*c2;
				bz1   = bz0*c1 + (bx0*By0 - by0*Bx0)*s1 + Bz0*bdotv*c2;

				ct0         = ct0 + (bz0/beta0 - 1)*cdt;

				x0  = x1;
				y0  = y1;
				z0  = z1;

				bx0 = bx1;
				by0 = by1;
				bz0 = bz1;
			}

			px1 = bx0*gamma/b0g0;
			py1 = by0*gamma/b0g0;
			
			beamPS[n*6]     = x0;
			beamPS[n*6 + 1] = px1;
			beamPS[n*6 + 2] = y0;
			beamPS[n*6 + 3] = py1;
			beamPS[n*6 + 4] = ct0;

			if(beamParams->g != 0) {
				theta = beamSpins[n*2];
				phi   = beamSpins[n*2 + 1];

				RotateSpins(Bx0, By0, Bz0, ds, px1, py1, dp0, &theta, &phi);

				beamSpins[n*2] = theta;
				beamSpins[n*2 + 1] = phi;
			}

			distance[n*2] += ds;

			if( (ax>0) & (ay>0) )
				if( (x0/ax)*(x0/ax) + (y0/ay)*(y0/ay) > 1)
					distance[n*2 + 1] = 0;
		}
	}

	beamParams->globaltime += fieldmap.length / (beta0*SpeedOfLight);

}
