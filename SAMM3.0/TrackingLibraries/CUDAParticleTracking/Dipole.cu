#include <cuda_runtime.h>
#include <math.h>

#include "CUDAParticleTracking.h"
#include "TrackCUDA.h"

#define _nparams 11
#define _length 0
#define _curvature 1
#define _field 2
#define _gradient 3
#define _e1 4
#define _e2 5
#define _hgap 6
#define _fint1 7
#define _fint2 8
#define _apertureX 9
#define _apertureY 10

extern "C" BeamParameters_t* beamParams;

extern Beam* d_BeamParameters;
extern Particle* d_ParticleArray;
extern int blocksPerGrid;

__global__ void DipoleMap(Component* dipole, Beam* beam, Particle* ptcle)
{
	int i = blockDim.x * blockIdx.x + threadIdx.x;

	const double brho   = beam[_rigidity];
	const double beta0  = beam[_beta];

	const double ds     = dipole[_length];
	const double h      = dipole[_curvature];
	const double e1     = dipole[_e1];
	const double e2     = dipole[_e2];
	const double hgap   = dipole[_hgap];
	const double k0     = dipole[_field] / brho;
	const double fint1  = dipole[_fint1];
	const double fint2  = dipole[_fint2];
	const double k1     = dipole[_gradient] / brho;
	const double ax     = dipole[_apertureX];
	const double ay     = dipole[_apertureY];

	const double sine1  = sin(e1);
	const double phi1   = 2*fint1*hgap*k0*(1+sine1*sine1)/cos(e1);
	const double sine2  = sin(e2);
	const double phi2   = 2*fint2*hgap*k0*(1+sine2*sine2)/cos(e2);
	
	const double r101   = k0*tan(e1);
    const double r321   =-k0*tan(e1 - phi1);

    const double r102   = k0*tan(e2);
    const double r322   =-k0*tan(e2 - phi2);
    
	
    if(ptcle[_f(i)])
	{
		double x0  = ptcle[_x(i)];
		double px0 = ptcle[_px(i)];
		double y0  = ptcle[_y(i)];
		double py0 = ptcle[_py(i)];
		double ct0 = ptcle[_ct(i)];
		double dp0 = ptcle[_dp(i)];
		
		double d1  = sqrt(1 + 2*dp0/beta0 + dp0*dp0);
	                
		double a1  = h - k0/d1;
	    
		double wx  = sqrt((h*k0 + k1)/d1);
		double xc  = cos(wx*ds);
		double xs  = sin(wx*ds)/wx;
		double xs2 = sin(2*wx*ds)/wx;
	    
		double wy  = sqrt(k1/d1);
		double yc  = cosh(wy*ds);
		double ys  = ds;
		double ys2 = 2*ds;

		if(wy!=0)
		{
			ys  = sinh(wy*ds)/wy;
			ys2 = sinh(2*wy*ds)/wy;
		}

		// First, apply a map for the entrance fringe field
		px0 = px0 + r101*x0/d1/d1;
		py0 = py0 + r321*y0/d1/d1;

		// Then, apply a map for the body of the dipole
		double x1 = x0*xc + px0*xs/d1 + a1*(1-xc)/wx/wx;
		double y1 = y0*yc + py0*ys/d1;

		double px1 =-d1*wx*wx*x0*xs + px0*xc + a1*xs*d1 ;
		double py1 = d1*wy*wy*y0*ys + py0*yc;

		double d0  = 1/beta0 + dp0;

		double c0  = (1/beta0 - d0/d1)*ds -
		             d0*a1*(h*(ds-xs) +
		                    a1*(2*ds-xs2)/8)/wx/wx/d1;
		
		double c1  =-d0*(h*xs - 
		                 a1*(2*ds-xs2)/4)/d1;
		
		double c2  =-d0*(h*(1-xc)/wx/wx +
		                 a1*xs*xs/2)/d1/d1;

		double c11 =-d0*wx*wx*(2*ds - xs2)/d1/8;
		double c12 = d0*wx*wx*xs*xs/d1/d1/2;
		double c22 =-d0*(2*ds + xs2)/d1/d1/d1/8;
		
		double c33 =-d0*wy*wy*(2*ds - ys2)/d1/8;
		double c34 =-d0*wy*wy*ys*ys/d1/d1/2;
		double c44 =-d0*(2*ds + ys2)/d1/d1/d1/8;
	    
		double ct1 = ct0 + c0 +
					c1*x0 + c2*px0 +
					c11*x0*x0 + c12*x0*px0 + c22*px0*px0 +
					c33*y0*y0 + c34*y0*py0 + c44*py0*py0;

		// Finally, apply a map for the exit fringe field
		px1 = px1 + r102*x1/d1/d1;
		py1 = py1 + r322*y1/d1/d1;

		ptcle[_x(i)]  = x1;
		ptcle[_px(i)] = px1;
		ptcle[_y(i)]  = y1;
		ptcle[_py(i)] = py1;
		ptcle[_ct(i)] = ct1;

		ptcle[_s(i)] += ds;
		
		if( (ax>0) & (ay>0) )
			if((x1/ax)*(x1/ax) + (y1/ay)*(y1/ay) > 1)
				ptcle[_f(i)] = 0;

	}

	if(i==0)
		beam[_globaltime] += ds / (beta0 * SpeedOfLight);

}

__global__ void DipoleSpinMap(Component* dipole, Beam* beam, Particle* ptcle)
{
	int i = blockDim.x * blockIdx.x + threadIdx.x;
	
	const double ds        = dipole[_length];
	const double field     = dipole[_field];
	const double gradient  = dipole[_gradient];
	const double bendangle = dipole[_length] * dipole[_curvature];
	
	const double gammaRef  = beam[_gamma];
	const double emratio   = beam[_charge] / beam[_mass];
	const double bigG      = (beam[_g] - 2)/2;
		
	if(ptcle[_f(i)])
	{
		double betaRef;
		
		double polsnx0;
		double polsny0;
		double polsnz0;

		double gamma;
		double beta;

		double ptot;
		double pz2;
		double pz;

		double pdotb;
		
		double bx;
		double by;
		double bz;

		double bParx;
		double bPary;
		double bParz;

		double bPerpx;
		double bPerpy;
		double bPerpz;

		double emdg;

		double omegax;
		double omegay;
		double omegaz;

		double omega;

		double pdotomega;

		double coswt;
		double sinwt;
		double c1;

		double polsnx1;
		double polsny1;
		double polsnz1;
		
		double px;
		double py;
		double dp;
		
		double theta;
		double phi;
		
		px       = ptcle[_px(i)];
		py       = ptcle[_py(i)];
		dp       = ptcle[_dp(i)];
		
		theta    = ptcle[_theta(i)];
		phi      = ptcle[_phi(i)];
		
		bx       = gradient*(ptcle[_y(i)]);
		by       = field + gradient*(ptcle[_x(i)]);
		bz       = 0;		

		polsnx0  = sin(theta)*cos(phi);
		polsny0  = sin(theta)*sin(phi);
		polsnz0  = cos(theta);

		betaRef  = sqrt(1 - 1/gammaRef/gammaRef);
		gamma    = (dp*betaRef + 1)*gammaRef;
		beta     = sqrt(1 - 1/gamma/gamma);
		ptot     = beta*gamma/betaRef/gammaRef;
		pz2      = ptot*ptot - px*px - py*py;
		if(pz2<0) pz2 = 0;
		pz       = sqrt(pz2);

		pdotb    = (px*bx + py*by + pz*bz)/ptot/ptot;

		bParx    = pdotb*px;
		bPary    = pdotb*py;
		bParz    = pdotb*pz;

		bPerpx   = bx - bParx;
		bPerpy   = by - bPary;
		bPerpz   = bz - bParz;
		
		emdg     = emratio / gamma;

		omegax   = -emdg*((1 + bigG*gamma)*bPerpx + (1 + bigG)*bParx);
		omegay   = -emdg*((1 + bigG*gamma)*bPerpy + (1 + bigG)*bPary);
		omegaz   = -emdg*((1 + bigG*gamma)*bPerpz + (1 + bigG)*bParz);

		omega    = sqrt(omegax*omegax + omegay*omegay + omegaz*omegaz);

		if(omega!=0) {
			pdotomega = polsnx0*omegax + polsny0*omegay + polsnz0*omegaz;

			coswt    = cos(omega*ds/SpeedOfLight);
			sinwt    = sin(omega*ds/SpeedOfLight);

			c1       = pdotomega*(1-coswt)/omega/omega;

			polsnx1  =  polsnx0*coswt +       c1*omegax +
					   (polsnz0*omegay - polsny0*omegaz)*sinwt/omega;

			polsny1  =  polsny0*coswt +       c1*omegay +
					   (polsnx0*omegaz - polsnz0*omegax)*sinwt/omega;

			polsnz1  =  polsnz0*coswt +       c1*omegaz +
					   (polsny0*omegax - polsnx0*omegay)*sinwt/omega;
		} else {
			polsnx1 = polsnx0;
			polsny1 = polsny0;
			polsnz1 = polsnz0;
		}
		
		theta   = acos(polsny1);
		phi     = atan2(polsnx1,polsnz1) + bendangle;

		polsnx1 = sin(theta)*sin(phi);
		polsny1 = cos(theta);
		polsnz1 = sin(theta)*cos(phi);

		ptcle[_theta(i)] = acos(polsnz1);
		ptcle[_phi(i)]   = atan2(polsny1,polsnx1);
	}
}

extern "C" __host__ void TrackDipoleCUDA(DipoleParameters_t dipole)
{
	size_t cptParamsListSize = _nparams*sizeof(Component);

	Component* d_ComponentParameters;

	Component parameters[_nparams];
	parameters[_length]    = dipole.length;
	parameters[_curvature] = dipole.curvature;
	parameters[_field]     = dipole.field;
	parameters[_gradient]  = dipole.gradient;
	parameters[_e1]        = dipole.e1;
	parameters[_e2]        = dipole.e2;
	parameters[_hgap]      = dipole.hgap;
	parameters[_fint1]     = dipole.fint1;
	parameters[_fint2]     = dipole.fint2;
	parameters[_apertureX] = dipole.apertureX;
	parameters[_apertureY] = dipole.apertureY;

	cudaMalloc((void**)&d_ComponentParameters, cptParamsListSize);
	cudaMemcpy(d_ComponentParameters, parameters, cptParamsListSize, cudaMemcpyHostToDevice);

	DipoleMap<<<blocksPerGrid, threadsPerBlock>>>(d_ComponentParameters, d_BeamParameters, d_ParticleArray);

	cudaFree(d_ComponentParameters);


	if(beamParams->g != 0)
	{
		cptParamsListSize = 4*sizeof(Component);

		Component diptclepinParams[4];
		diptclepinParams[_length]    = dipole.length;
		diptclepinParams[_field]     = dipole.field;
		diptclepinParams[_gradient]  = dipole.gradient;
		diptclepinParams[_curvature] = dipole.curvature;

		cudaMalloc((void**)&d_ComponentParameters, cptParamsListSize);
		cudaMemcpy(d_ComponentParameters, diptclepinParams, cptParamsListSize, cudaMemcpyHostToDevice);

		DipoleSpinMap<<<blocksPerGrid, threadsPerBlock>>>(d_ComponentParameters, d_BeamParameters, d_ParticleArray);

		cudaFree(d_ComponentParameters);
	}
}
