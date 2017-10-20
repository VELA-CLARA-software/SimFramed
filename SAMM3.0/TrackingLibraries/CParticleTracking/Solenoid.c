
#include <math.h>

#include "CParticleTracking.h"
#include "GlobalVariables.h"

void TrackSolenoid(SolenoidParameters_t solenoid)
{
	const double ds     = solenoid.length;
	const double b0     = solenoid.field;
    const double g      = solenoid.taper;
            
    const double brho   = beamParams->rigidity;
	const double b0g0   = brho * beamParams->charge / beamParams->mass / SpeedOfLight;
	const double gamma0 = sqrt(1+b0g0*b0g0);
	const double beta0  = b0g0/gamma0;
			
	const double ax = solenoid.apertureX;
	const double ay = solenoid.apertureY;

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


			const double gds  = g*ds;
			const double b1   = b0 / (1 + gds);
            
            const double helmA = sqrt(b0/b1);
            const double helmB = 2*brho*(1+dp0)/sqrt(b0*b1);
            const double helmF = -1 / helmB;
            const double helmG =  1 / helmA;
            
            double w = b0*ds/2/brho/(1+dp0);

			double cw2;
			double s2w;
			double sw2;

			double x1;
			double px1;
			double y1;
			double py1;

			double theta;
			double phi;

			double bx;
			double by;
			double bz;

			double d1  = sqrt(1 - px0*px0 - py0*py0
				                + 2*dp0/beta0 + dp0*dp0);

			ct0 += ds*(1 - (1 + beta0*dp0)/d1)/beta0/2;

            if(gds!=0)
                w = w*log(1+gds)/gds;
            
            cw2 = cos(w)*cos(w);
            s2w = sin(2*w);
            sw2 = sin(w)*sin(w);
            
            x1  = helmA*cw2*x0   + helmB*s2w*px0/2 + helmA*s2w*y0/2 + helmB*sw2*py0;
            px1 = helmF*s2w*x0/2 + helmG*cw2*px0   + helmF*sw2*y0   + helmG*s2w*py0/2;
            y1  =-helmA*s2w*x0/2 - helmB*sw2*px0   + helmA*cw2*y0   + helmB*s2w*py0/2;
            py1 =-helmF*sw2*x0   - helmG*s2w*px0/2 + helmF*s2w*y0/2 + helmG*cw2*py0;


			d1  = sqrt(1 - px1*px1 - py1*py1
			             + 2*dp0/beta0 + dp0*dp0);

			ct0 += ds*(1 - (1 + beta0*dp0)/d1)/beta0/2;

			beamPS[n*6]     = x1;
			beamPS[n*6 + 1] = px1;
			beamPS[n*6 + 2] = y1;
			beamPS[n*6 + 3] = py1;
			beamPS[n*6 + 4] = ct0;

			if(beamParams->g != 0) {
				theta = beamSpins[n*2];
				phi   = beamSpins[n*2 + 1];

				bx    = x1*b1*g/2;
				by    = y1*b1*g/2;
				bz    = b0;

				if(gds != 0)
					bz = bz*log(1+gds)/gds;

				RotateSpins(bx, by, bz, ds, px1, py1, dp0, &theta, &phi);

				beamSpins[n*2] = theta;
				beamSpins[n*2 + 1] = phi;
			}

			distance[n*2] += ds;

			if( (ax>0) & (ay>0) )
				if( (x1/ax)*(x1/ax) + (y1/ay)*(y1/ay) > 1)
					distance[n*2 + 1] = 0;
		}
	}

	beamParams->globaltime += solenoid.length / (beta0*SpeedOfLight);

}
