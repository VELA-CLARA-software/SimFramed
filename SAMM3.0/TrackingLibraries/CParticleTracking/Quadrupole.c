
#include <math.h>
#include <stdlib.h>

#include "CParticleTracking.h"
#include "AppendComponent.h"
#include "GlobalVariables.h"

#define _length 0
#define _gradient 1
#define _apertureX 2
#define _apertureY 3


void TrackQuadrupole_(double* parameters)
{
	const double ds = parameters[_length];
	const double k1 = parameters[_gradient] / beamParams->rigidity;

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

			double d1  = sqrt(1 + 2*dp0/beta0 + dp0*dp0);

			double w   = sqrt(fabs(k1)/d1);
	        
			double xs;
			double xc;
			double ys;
			double yc;
			double xs2;
			double ys2;

			double x1;
			double px1;
			double y1;
			double py1;

			double d0;
			double d2;

			double c0;
			double c11;
			double c12;
			double c22;
			double c33;
			double c34;
			double c44;

			double ct1;

			double theta;
			double phi;

			double bx;
			double by;
	        
			if(k1>=0)
			{
			   xs  = sin(w*ds);
			   xc  = cos(w*ds);
			   ys  = sinh(w*ds);
			   yc  = cosh(w*ds);
			   xs2 = sin(2*w*ds);
			   ys2 = sinh(2*w*ds);
			}
			else
			{
			   xs  = sinh(w*ds);
			   xc  = cosh(w*ds);
			   ys  = sin(w*ds);
			   yc  = cos(w*ds);
			   xs2 = sinh(2*w*ds);
			   ys2 = sin(2*w*ds);
			}
	        
			x1  =  x0*xc      + px0*xs*w/fabs(k1);
			px1 = -k1*x0*xs/w + px0*xc;
			y1  =  y0*yc      + py0*ys*w/fabs(k1);
			py1 =  k1*y0*ys/w + py0*yc;

			d0  = 1/beta0 + dp0;
			d2  =-d0/d1/d1/d1/2;

			c0  = (1/beta0 - d0/d1)*ds;
	        c11 = k1*k1*d2*(xs2/w - 2*ds)/w/w/4;
			c12 =-k1*d2*xs*xs/w/w;
			c22 = d2*(xs2/w + 2*ds)/4;
			c33 = k1*k1*d2*(ys2/w - 2*ds)/w/w/4;
			c34 = k1*d2*ys*ys/w/w;
			c44 = d2*(ys2/w + 2*ds)/4;
	        
			ct1 = ct0 + c0
				      + c11* x0* x0
				      + c12* x0*px0
				      + c22*px0*px0
				      + c33* y0* y0
				      + c34* y0*py0
				      + c44*py0*py0;

			beamPS[n*6]     = x1;
			beamPS[n*6 + 1] = px1;
			beamPS[n*6 + 2] = y1;
			beamPS[n*6 + 3] = py1;
			beamPS[n*6 + 4] = ct1;

			if(beamParams->g != 0) {
				theta = beamSpins[n*2];
				phi   = beamSpins[n*2 + 1];

				bx    = parameters[_gradient] * y1;
				by    = parameters[_gradient] * x1;

				RotateSpins(bx, by, 0, ds, px1, py1, dp0, &theta, &phi);

				beamSpins[n*2] = theta;
				beamSpins[n*2 + 1] = phi;
			}

			distance[n*2] += ds;

			if( (ax>0) & (ay>0) )
				if( (x1/ax)*(x1/ax) + (y1/ay)*(y1/ay) > 1)
					distance[n*2 + 1] = 0;
		}
	}

	beamParams->globaltime += parameters[_length] / (beta0*SpeedOfLight);

}


void TrackQuadrupole(QuadrupoleParameters_t quadrupole)
{
	double parameters[4];

	parameters[_length]    = quadrupole.length;
	parameters[_gradient]  = quadrupole.gradient;
	parameters[_apertureX] = quadrupole.apertureX;
	parameters[_apertureY] = quadrupole.apertureY;

	TrackQuadrupole_(parameters);
}


void CopyQuadrupole(QuadrupoleParameters_t quadrupole)
{
	double* parameters = (double*)calloc(4,sizeof(double));

	parameters[_length]    = quadrupole.length;
	parameters[_gradient]  = quadrupole.gradient;
	parameters[_apertureX] = quadrupole.apertureX;
	parameters[_apertureY] = quadrupole.apertureY;

	AppendComponent(TrackQuadrupole_, parameters);
}