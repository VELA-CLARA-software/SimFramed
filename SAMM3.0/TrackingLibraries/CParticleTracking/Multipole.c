
#include <math.h>
#include <stdlib.h>

#include "CParticleTracking.h"
#include "AppendComponent.h"
#include "GlobalVariables.h"

#define _apertureX 0
#define _apertureY 1
#define _bendangle 2
#define _curvature 3
#define _ncpts 4
#define _findex 5


// int factorial(int n)
// {
//	return (n==1 || n==0) ? 1 : factorial(n - 1) * n;
// }


void TrackMultipole_(double* parameters)
{
	const double ax        = parameters[_apertureX];
	const double ay        = parameters[_apertureY];

	const double bendangle = parameters[_bendangle];
	const double curvature = parameters[_curvature];

	const int    ncpts     = (int)parameters[_ncpts];

	double*      findex    = &parameters[_findex];
	double*      anL       = &parameters[_findex + ncpts];
	double*      bnL       = &parameters[_findex + ncpts*2];

	const double brho      = beamParams->rigidity;

	int          n;


	for(n=0; n<beamParams->nparticles; n++)
	{
		if(distance[n*2 + 1])
		{

			double x0   = beamPS[n*6];
			double px0  = beamPS[n*6 + 1];
			double y0   = beamPS[n*6 + 2];
			double py0  = beamPS[n*6 + 3];
			double ct0  = beamPS[n*6 + 4];
			double dp0  = beamPS[n*6 + 5];

			double BxL  = 0;
			double ByL  = 0;

			double psi  = atan2(y0,x0);
			double r    = sqrt(x0*x0 + y0*y0);

			int    indx = 0;
			int    ni   = 1;

			double pr   = 1;

			for(indx=0; indx<ncpts; indx++)
			{
				int nn    = (int) findex[indx];

				double cs = cos(nn*psi);
				double sn = sin(nn*psi);

				//double pr = pow(r,nn) / factorial(nn);
				for( ; ni<=nn; ni++)
					pr = pr * r / ni;

				ByL      += pr*(bnL[indx]*cs - anL[indx]*sn);
				BxL      += pr*(bnL[indx]*sn + anL[indx]*cs);
			}

			px0 += bendangle*(1 + dp0 - curvature*x0) - ByL/brho;
			py0 += BxL/brho;
			ct0 -= bendangle*x0;

			beamPS[n*6 + 1] = px0;
			beamPS[n*6 + 3] = py0;
			beamPS[n*6 + 4] = ct0;


			if(beamParams->g != 0) {
				double theta = beamSpins[n*2];
				double phi   = beamSpins[n*2 + 1];
				double L     = 1e-9;
				double BzL   = 0;

				RotateSpins(BxL/L, ByL/L, BzL/L, L, px0, py0, dp0, &theta, &phi);

				beamSpins[n*2] = theta;
				beamSpins[n*2 + 1] = phi;
			}


			if( (ax>0) & (ay>0) )
				if( (x0/ax)*(x0/ax) + (y0/ay)*(y0/ay) > 1)
					distance[n*2 + 1] = 0;

		}
	}

}


void TrackMultipole(MultipoleParameters_t multipole)
{
	int ncpts  = (int) multipole.ncpts;

	int n;

	double* parameters = (double*)calloc(5+3*ncpts,sizeof(double));

	parameters[_apertureX] = multipole.apertureX;
	parameters[_apertureY] = multipole.apertureY;

	parameters[_bendangle] = multipole.angle;
	parameters[_curvature] = multipole.curvature;

	parameters[_ncpts]     = (double)ncpts;

	for(n=0; n<ncpts; n++)
	{
		parameters[_findex + n]            = (double)multipole.fieldindex[n];
		parameters[_findex + n + ncpts]    = multipole.anL[n];
		parameters[_findex + n + ncpts*2]  = multipole.bnL[n];
	}

	TrackMultipole_(parameters);

	free(parameters);
}


void CopyMultipole(MultipoleParameters_t multipole)
{
	int ncpts  = (int) multipole.ncpts;

	int n;

	double* parameters = (double*)calloc(5+3*ncpts,sizeof(double));

	parameters[_apertureX] = multipole.apertureX;
	parameters[_apertureY] = multipole.apertureY;

	parameters[_bendangle] = multipole.angle;
	parameters[_curvature] = multipole.curvature;

	parameters[_ncpts]     = (double)ncpts;

	for(n=0; n<ncpts; n++)
	{
		parameters[_findex + n]            = (double)multipole.fieldindex[n];
		parameters[_findex + n + ncpts]    = multipole.anL[n];
		parameters[_findex + n + ncpts*2]  = multipole.bnL[n];
	}

	AppendComponent(TrackMultipole_, parameters);
}
