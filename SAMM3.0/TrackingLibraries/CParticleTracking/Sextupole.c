
#include <math.h>

#include "CParticleTracking.h"
#include "GlobalVariables.h"

void TrackSextupole(SextupoleParameters_t sextupole)
{
	const double ds = sextupole.length;
	const double k2 = sextupole.gradient / beamParams->rigidity;

	const double ax = sextupole.apertureX;
	const double ay = sextupole.apertureY;

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

			double theta;
			double phi;

			double bx;
			double by;

			// First, apply a drift through ds/2
			double d1  = sqrt(1 - px0*px0 - py0*py0
				                + 2*dp0/beta0 + dp0*dp0);

			x0  += ds*px0/d1/2;
			y0  += ds*py0/d1/2;
			ct0 += ds*(1 - (1 + beta0*dp0)/d1)/beta0/2;

			// Then, apply a sextupole kick
            px0 += -(x0*x0 - y0*y0)*k2*ds/2;
            py0 += x0*y0*k2*ds;

			// Next, apply another drift through ds/2
			d1  = sqrt(1 - px0*px0 - py0*py0
			             + 2*dp0/beta0 + dp0*dp0);

			x0  += ds*px0/d1/2;
			y0  += ds*py0/d1/2;
			ct0 += ds*(1 - (1 + beta0*dp0)/d1)/beta0/2;

			beamPS[n*6]     = x0;
			beamPS[n*6 + 1] = px0;
			beamPS[n*6 + 2] = y0;
			beamPS[n*6 + 3] = py0;
			beamPS[n*6 + 4] = ct0;

			if(beamParams->g != 0) {
				theta = beamSpins[n*2];
				phi   = beamSpins[n*2 + 1];

				bx    = sextupole.gradient * (x0*x0 - y0*y0)/2;
				by    = sextupole.gradient * x0 * y0;

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

	beamParams->globaltime += sextupole.length / (beta0*SpeedOfLight);

}
