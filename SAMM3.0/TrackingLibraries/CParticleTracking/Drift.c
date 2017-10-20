
#include <math.h>
#include <stdlib.h>

#include "CParticleTracking.h"
#include "AppendComponent.h"
#include "GlobalVariables.h"

#define _length 0
#define _apertureX 1
#define _apertureY 2


void TrackDrift_(double* parameters)
{
	const double ds = parameters[_length];
	const double ax = parameters[_apertureX];
	const double ay = parameters[_apertureY];

	const double b0g0   = beamParams->rigidity * beamParams->charge / beamParams->mass / SpeedOfLight;
	const double gamma0 = sqrt(1+b0g0*b0g0);
	const double beta0  = b0g0/gamma0;

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

			double d1  = sqrt(1 - px0*px0 - py0*py0
				                + 2*dp0/beta0 + dp0*dp0);

			x0  += ds*px0/d1;
			y0  += ds*py0/d1;
			ct0 += ds*(1 - (1 + beta0*dp0)/d1)/beta0;

			beamPS[n*6]     = x0;
			beamPS[n*6 + 2] = y0;
			beamPS[n*6 + 4] = ct0;

			distance[n*2] += ds;

			if( (ax>0) & (ay>0) )
				if( (x0/ax)*(x0/ax) + (y0/ay)*(y0/ay) > 1)
					distance[n*2 + 1] = 0;
		}
	}

	beamParams->globaltime += ds / (beta0*SpeedOfLight);
}


void TrackDrift(DriftParameters_t drift)
{
	double parameters[3];

	parameters[_length]    = drift.length;
	parameters[_apertureX] = drift.apertureX;
	parameters[_apertureY] = drift.apertureY;

	TrackDrift_(parameters);
}


void CopyDrift(DriftParameters_t drift)
{
	double* parameters = (double*)calloc(3,sizeof(double));

	parameters[_length]    = drift.length;
	parameters[_apertureX] = drift.apertureX;
	parameters[_apertureY] = drift.apertureY;

	AppendComponent(TrackDrift_, parameters);
}