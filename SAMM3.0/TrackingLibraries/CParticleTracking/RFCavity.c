
#include <math.h>
#include <stdlib.h>

#include "CParticleTracking.h"
#include "AppendComponent.h"
#include "GlobalVariables.h"

#define pi (double)3.1415926536

#define _length 0
#define _voltage 1
#define _frequency 2
#define _phase 3
#define _apertureX 4
#define _apertureY 5
#define _masteroscillatorfrequency 6


void TrackRFCavity_(double* parameters)
{
	const double ds  = parameters[_length];
	const double v0  = parameters[_voltage] / beamParams->rigidity / SpeedOfLight;
	const double f   = parameters[_frequency];
	const double phi = parameters[_phase];
	const double ax  = parameters[_apertureX];
	const double ay  = parameters[_apertureY];
	const double mo  = parameters[_masteroscillatorfrequency];

	const double b0g0   = beamParams->rigidity * beamParams->charge / beamParams->mass / SpeedOfLight;
	const double gamma0 = sqrt(1+b0g0*b0g0);
	const double beta0  = b0g0/gamma0;

	int n;

	int p = (int)(beamParams->globaltime * mo);
	beamParams->globaltime -= (double)p / mo;

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

			// First, apply a drift through ds/2
			double d1  = sqrt(1 - px0*px0 - py0*py0
				                + 2*dp0/beta0 + dp0*dp0);

			double t;

			x0  += ds*px0/d1/2;
			y0  += ds*py0/d1/2;
			ct0 += ds*(1 - (1 + beta0*dp0)/d1)/beta0/2;

			// Then, apply an RF kick
			t = beamParams->globaltime - ct0 / (beta0*SpeedOfLight);
			dp0 += v0*sin(2*pi*f*t + phi);

			// Next, apply another drift through ds/2
			d1  = sqrt(1 - px0*px0 - py0*py0
				         + 2*dp0/beta0 + dp0*dp0);

			x0  += ds*px0/d1/2;
			y0  += ds*py0/d1/2;
			ct0 += ds*(1 - (1 + beta0*dp0)/d1)/beta0/2;

			beamPS[n*6]     = x0;
			beamPS[n*6 + 2] = y0;
			beamPS[n*6 + 4] = ct0;
			beamPS[n*6 + 5] = dp0;

			// Finally, collimate
			distance[n*2] += ds;

			if( (ax>0) & (ay>0) )
				if( (x0/ax)*(x0/ax) + (y0/ay)*(y0/ay) > 1)
					distance[n*2 + 1] = 0;
		}
	}

	beamParams->globaltime += ds / (beta0*SpeedOfLight);
}


void TrackRFCavity(RFCavityParameters_t rfcavity)
{
	double parameters[7];

	parameters[_length]                     = rfcavity.length;
	parameters[_voltage]                    = rfcavity.voltage;
	parameters[_frequency]                  = rfcavity.frequency;
	parameters[_phase]                      = rfcavity.phase;
	parameters[_apertureX]                  = rfcavity.apertureX;
	parameters[_apertureY]                  = rfcavity.apertureY;
	parameters[_masteroscillatorfrequency]  = rfcavity.masteroscillatorfrequency;

	TrackRFCavity_(parameters);
}


void CopyRFCavity(RFCavityParameters_t rfcavity)
{
	double* parameters = (double*)calloc(7,sizeof(double));

	parameters[_length]                     = rfcavity.length;
	parameters[_voltage]                    = rfcavity.voltage;
	parameters[_frequency]                  = rfcavity.frequency;
	parameters[_phase]                      = rfcavity.phase;
	parameters[_apertureX]                  = rfcavity.apertureX;
	parameters[_apertureY]                  = rfcavity.apertureY;
	parameters[_masteroscillatorfrequency]  = rfcavity.masteroscillatorfrequency;

	AppendComponent(TrackRFCavity_, parameters);
}
