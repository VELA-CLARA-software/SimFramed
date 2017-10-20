#include <stdlib.h>
#include <math.h>

#include "CParticleTracking.h"
#include "AppendComponent.h"
#include "GlobalVariables.h"

#define _length 0
#define _apertureX 1
#define _apertureY 2
#define _ncoeffs 3
#define _coeffs 4

void TrackTaylorMap_(double* parameters)
{

	const double ax     = parameters[_apertureX];
	const double ay     = parameters[_apertureY];

	const double b0g0   = beamParams->rigidity * beamParams->charge / beamParams->mass / SpeedOfLight;
	const double beta0  = b0g0/sqrt(1+b0g0*b0g0);

	int n;

	for(n=0; n<beamParams->nparticles; n++)
	{
		double x0  = beamPS[n*6];
		double y0  = beamPS[n*6 + 2];

		if( (ax>0) & (ay>0) )
			if( (x0/ax)*(x0/ax) + (y0/ay)*(y0/ay) > 1)
				distance[n*2 + 1] = 0;

		if(distance[n*2 + 1])
		{
			double ps0[6];

			int m, nc;

			for(m=0; m<6; m++)
			{
				ps0[m] = beamPS[6*n + m];
				beamPS[6*n + m] = 0;
			}

			for(nc=0; nc<parameters[_ncoeffs]; nc++)
			{
				int var  = (int) parameters[_coeffs + nc*8];
				double c =       parameters[_coeffs + nc*8 + 1];

				for(m=0; m<6; m++)
					c *=  pow(ps0[m], parameters[_coeffs + nc*8 + 2 + m]);

				beamPS[6*n + var - 1] += c;
			}

		}

	}

	beamParams->globaltime += parameters[_length] / (beta0*SpeedOfLight);

}


void TrackTaylorMap(TaylorMapParameters_t tmap)
{  
	int n;

	double* parameters = (double*) calloc(4 + tmap.ncoeffs*8, sizeof(double));
  
	parameters[_length]    = tmap.length;
	parameters[_apertureX] = tmap.apertureX;
	parameters[_apertureY] = tmap.apertureY;
	parameters[_ncoeffs]   = (double) tmap.ncoeffs;

	for (n=0; n<tmap.ncoeffs*8; n++)
		parameters[_coeffs + n]  = tmap.coeffs[n];

	TrackTaylorMap_(parameters);
  
	free(parameters);
}


void CopyTaylorMap(TaylorMapParameters_t tmap)
{
	int n;

	double* parameters = (double*) calloc(4 + tmap.ncoeffs*8, sizeof(double));
  
	parameters[_length]    = tmap.length;
	parameters[_apertureX] = tmap.apertureX;
	parameters[_apertureY] = tmap.apertureY;
	parameters[_ncoeffs]   = (double) tmap.ncoeffs;

	for (n=0; n<tmap.ncoeffs*8; n++)
		parameters[_coeffs + n]  = tmap.coeffs[n];

	AppendComponent(TrackTaylorMap_, parameters);
}
