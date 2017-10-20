
#include <math.h>
#include <stdlib.h>

#include "CParticleTracking.h"
#include "AppendComponent.h"
#include "GlobalVariables.h"

#define _length 0
#define _fieldX 1
#define _fieldY 2
#define _apertureX 3
#define _apertureY 4


void TrackOrbitCorrector_(double* parameters)
{
	const double ds  = parameters[_length];
	
	const double bx  = parameters[_fieldX];
	const double by  = parameters[_fieldY];

	const double k0x = bx / beamParams->rigidity;
	const double k0y = by / beamParams->rigidity;

	const double ax  = parameters[_apertureX];
	const double ay  = parameters[_apertureY];

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

			double x1;
			double px1;
			double y1;
			double py1;
			double ct1;

			double f1;
			double c0;

			double theta;
			double phi;

			double d1  = sqrt(1 + 2*dp0/beta0 + dp0*dp0);

			x1  = x0 + ds*px0/d1 - ds*ds*k0y/d1/2;
			px1 =         px0    - ds*k0y;

			y1  = y0 + ds*py0/d1 + ds*ds*k0x/d1/2;
			py1 =         py0    + ds*k0x;

			f1  = ds*(1/beta0 + dp0)/d1/d1/d1/2;

			c0  = ct0 + ds/beta0 - ds*(1/beta0 + dp0)/d1 -
				  ds*ds*f1*(k0x*k0x + k0y*k0y)/3;

			ct1 = c0 + ds*f1*(k0y*px0 - k0x*py0) - f1*(px0*px0 + py0*py0);

			beamPS[n*6]     = x1;
			beamPS[n*6 + 1] = px1;
			beamPS[n*6 + 2] = y1;
			beamPS[n*6 + 3] = py1;
			beamPS[n*6 + 4] = ct1;

			if(beamParams->g != 0) {
				theta = beamSpins[n*2];
				phi   = beamSpins[n*2 + 1];

				RotateSpins(bx, by, 0, ds, px0, py0, dp0, &theta, &phi);

				beamSpins[n*2] = theta;
				beamSpins[n*2 + 1] = phi;
			}

			// Finally, collimate
			distance[n*2] += ds;

			if( (ax>0) & (ay>0) )
				if( (x0/ax)*(x0/ax) + (y0/ay)*(y0/ay) > 1)
					distance[n*2 + 1] = 0;
		}
	}

	beamParams->globaltime += ds / (beta0*SpeedOfLight);

}


void TrackOrbitCorrector(OrbitCorrectorParameters_t orbitcorrector)
{
	double parameters[5];

	parameters[_length]    = orbitcorrector.length;
	parameters[_fieldX]    = orbitcorrector.fieldX;
	parameters[_fieldY]    = orbitcorrector.fieldY;
	parameters[_apertureX] = orbitcorrector.apertureX;
	parameters[_apertureY] = orbitcorrector.apertureY;

	TrackOrbitCorrector_(parameters);
}


void CopyOrbitCorrector(OrbitCorrectorParameters_t orbitcorrector)
{
	double* parameters = (double*)calloc(5,sizeof(double));

	parameters[_length]    = orbitcorrector.length;
	parameters[_fieldX]    = orbitcorrector.fieldX;
	parameters[_fieldY]    = orbitcorrector.fieldY;
	parameters[_apertureX] = orbitcorrector.apertureX;
	parameters[_apertureY] = orbitcorrector.apertureY;

	AppendComponent(TrackOrbitCorrector_, parameters);
}