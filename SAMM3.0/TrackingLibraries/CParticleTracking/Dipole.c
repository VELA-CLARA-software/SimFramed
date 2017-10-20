
#include <math.h>

#include "CParticleTracking.h"
#include "GlobalVariables.h"

void TrackDipole(DipoleParameters_t dipole)
{
	const double ds     = dipole.length;
	const double k0     = dipole.field / beamParams->rigidity;
    const double h      = dipole.curvature;
    const double k1     = dipole.gradient / beamParams->rigidity;

	const double sine1  = sin(dipole.e1);
	const double phi1   = 2*dipole.fint1*dipole.hgap*k0*(1+sine1*sine1)/cos(dipole.e1);
	const double sine2  = sin(dipole.e2);
	const double phi2   = 2*dipole.fint2*dipole.hgap*k0*(1+sine2*sine2)/cos(dipole.e2);

    const double r101   = k0*tan(dipole.e1);
    const double r321   =-k0*tan(dipole.e1 - phi1);
    
    const double r102   = k0*tan(dipole.e2);
    const double r322   =-k0*tan(dipole.e2 - phi2);

	const double ax     = dipole.apertureX;
	const double ay     = dipole.apertureY;

	const double b0g0   = beamParams->rigidity * beamParams->charge / beamParams->mass / SpeedOfLight;
	const double gamma0 = sqrt(1+b0g0*b0g0);
	const double beta0  = b0g0/gamma0;

	int n;
	for(n=0; n<beamParams->nparticles; n++)
	{
		if(distance[n*2 + 1])
		{
			double a1;
			double wx;
			double xc;
			double xs;
			double xs2;

			double wy;
			double yc;
			double ys;
			double ys2;

			double x1;
			double px1;
			double y1;
			double py1;
			double ct1;

			double d0;

			double c0;
			double c1;
			double c2;
			double c11;
			double c12;
			double c22;
			double c33;
			double c34;
			double c44;

			double x0  = beamPS[n*6];
			double px0 = beamPS[n*6 + 1];
			double y0  = beamPS[n*6 + 2];
			double py0 = beamPS[n*6 + 3];
			double ct0 = beamPS[n*6 + 4];
			double dp0 = beamPS[n*6 + 5];

			double d1  = sqrt(1 + 2*dp0/beta0 + dp0*dp0);

			double theta;
			double phi;

			double polsnx0;
			double polsny0;
			double polsnz0;
		                
			// First, apply a map for the entrance fringe field
			px0 = px0 + r101*x0/d1/d1;
			py0 = py0 + r321*y0/d1/d1;
		    
			// Then, apply a map for the body of the dipole
			a1  = h - k0/d1;
		    
			wx  = sqrt((h*k0 + k1)/d1);
			xc  = cos(wx*ds);
			xs  = sin(wx*ds)/wx;
			xs2 = sin(2*wx*ds)/wx;
		    
			wy  = sqrt(k1/d1);
			yc  = coshf(wy*ds);
			ys  = ds;
			ys2 = 2*ds;

			if(wy!=0)
			{
				ys  = sinhf(wy*ds)/wy;
				ys2 = sinhf(2*wy*ds)/wy;
			}
		    
			x1  =          x0*xc + px0*xs/d1 + a1*(1-xc)/wx/wx;
			px1 =-d1*wx*wx*x0*xs + px0*xc    + a1*xs*d1;
		    
			y1  =          y0*yc + py0*ys/d1;
			py1 = d1*wy*wy*y0*ys + py0*yc;

			d0  = 1/beta0 + dp0;

			c0  = (1/beta0 - d0/d1)*ds -
				  d0*a1*(h*(ds-xs) +
				  a1*(2*ds-xs2)/8)/wx/wx/d1;

			c1  =-d0*(h*xs - 
				  a1*(2*ds-xs2)/4)/d1;

			c2  =-d0*(h*(1-xc)/wx/wx +
				  a1*xs*xs/2)/d1/d1;

			c11 =-d0*wx*wx*(2*ds - xs2)/d1/8;
			c12 = d0*wx*wx*xs*xs/d1/d1/2;
			c22 =-d0*(2*ds + xs2)/d1/d1/d1/8;

			c33 =-d0*wy*wy*(2*ds - ys2)/d1/8;
			c34 =-d0*wy*wy*ys*ys/d1/d1/2;
			c44 =-d0*(2*ds + ys2)/d1/d1/d1/8;
		    
			ct1 = ct0 + c0 +
				  c1*x0 + c2*px0 +
				  c11*x0*x0 + c12*x0*px0 + c22*px0*px0 +
				  c33*y0*y0 + c34*y0*py0 + c44*py0*py0;
		    
			// Finally, apply a map for the exit fringe field
			px1 = px1 + r102*x1/d1/d1;
			py1 = py1 + r322*y1/d1/d1;

			beamPS[n*6]     = x1;
			beamPS[n*6 + 1] = px1;
			beamPS[n*6 + 2] = y1;
			beamPS[n*6 + 3] = py1;
			beamPS[n*6 + 4] = ct1;

			if(beamParams->g != 0) {
				theta = beamSpins[n*2];
				phi   = beamSpins[n*2 + 1];
				RotateSpins(0, dipole.field, 0, ds, px1, py1, dp0, &theta, &phi);

				polsnx0 = sin(theta)*cos(phi);
				polsny0 = sin(theta)*sin(phi);
				polsnz0 = cos(theta);

				theta   = acos(polsny0);
				phi     = atan2(polsnx0,polsnz0) + dipole.length * dipole.curvature;

				polsnx0 = sin(theta)*sin(phi);
				polsny0 = cos(theta);
				polsnz0 = sin(theta)*cos(phi);

				beamSpins[n*2] = acos(polsnz0);
				beamSpins[n*2 + 1] = atan2(polsny0,polsnx0);
			}

			distance[n*2] += ds;

			if( (ax>0) & (ay>0) )
				if( (x1/ax)*(x1/ax) + (y1/ay)*(y1/ay) > 1)
					distance[n*2 + 1] = 0;
		}
	}

	beamParams->globaltime += dipole.length / (beta0*SpeedOfLight);

}
