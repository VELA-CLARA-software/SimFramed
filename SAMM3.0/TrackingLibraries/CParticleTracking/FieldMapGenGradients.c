
#include <math.h>
#include <stdlib.h>

#include "CParticleTracking.h"
#include "AppendComponent.h"
#include "GlobalVariables.h"

#define SWAP(a,b) {temp=(a);(a)=(b);(b)=temp;}

#define _apertureX 0
#define _apertureY 1
#define _msize 2
#define _ssize 3
#define _lsize 4
#define _mvalues 5
#define _svalues (5 + msize)
#define _ctable (5 + msize + ssize)


void gaussj(double* a, double* b)
{
	int indxc[6];
	int indxr[6];
	int ipiv[6];

	int i;
	int icol;
	int irow;
	int j;
	int k;
	int l;
	int ll;

	double big;
	double dum;
	double pivinv;
	double temp;

	const int n = 6;

	for(j=0; j<n; j++)
		ipiv[j] = 0;

	for(i=0; i<n; i++)
	{
		big = 0.0;

		for(j=0; j<n; j++)
			if (ipiv[j] != 1)
				for(k=0; k<n; k++)
					if (ipiv[k] == 0)
						if (fabs(a[6*j + k]) >= big)
						{
							big = fabs(a[6*j + k]);
							irow = j;
							icol = k;
						}

		++(ipiv[icol]);

		if (irow != icol)
		{
			for (l=0; l<n; l++)
				SWAP(a[6*irow + l],a[6*icol + l]);

			SWAP(b[irow],b[icol]);
		}

		indxr[i] = irow;
		indxc[i] = icol;

		pivinv = 1.0/a[6*icol + icol];
		a[6*icol + icol] = 1.0;

		for (l=0; l<n; l++)
			a[6*icol + l] *= pivinv;
		b[icol] *= pivinv;

		for (ll=0; ll<n; ll++)
			if (ll != icol)
			{
				dum = a[6*ll + icol];
				a[6*ll + icol] = 0.0;

				for (l=0; l<n; l++)
					a[6*ll + l] -= a[6*icol + l]*dum;
				b[ll] -= b[icol]*dum;
			}
	}

	for (l=n-1; l>=0; l--)
		if (indxr[l] != indxc[l])
			for (k=0; k<n; k++)
				SWAP(a[6*k + indxr[l]],a[6*k + indxc[l]]);
}


void gradH(double* dH, double* avector, double* davector, double* ps0)
{
	const double ax     = avector[0];
	const double ay     = avector[1];

	const double axx    = davector[0];
	const double axy    = davector[1];
	const double ayx    = davector[2];
	const double ayy    = davector[3];
	const double asx    = davector[4];
	const double asy    = davector[5];

	const double qdmc   = beamParams->charge / beamParams->mass / SpeedOfLight;
	const double b0g0   = beamParams->rigidity * qdmc;
	const double gamma0 = sqrt(1+b0g0*b0g0);
	const double beta0  = b0g0/gamma0;

	const double px1    = ps0[1] - ax;
	const double py1    = ps0[3] - ay;
	const double dp1    = 1/beta0 + ps0[5];

	const double h1     = sqrt( dp1*dp1 - px1*px1 - py1*py1 - 1/b0g0/b0g0 );

	dH[0]               = -(px1*axx + py1*ayx)/h1 - asx;
	dH[1]               = px1/h1;

	dH[2]               = -(px1*axy + py1*ayy)/h1 - asy;
	dH[3]               = py1/h1;

	dH[4]               = 0;
	dH[5]               = 1/beta0 - dp1/h1;
}


void gradgradH(double* d2H, double* avector, double* davector, double* d2avector, double* ps0)
{
	const double ax     = avector[0];
	const double ay     = avector[1];

	const double axx    = davector[0];
	const double axy    = davector[1];
	const double ayx    = davector[2];
	const double ayy    = davector[3];
	const double asx    = davector[4];
	const double asy    = davector[5];

	const double axxx   = d2avector[0];
	const double axxy   = d2avector[1];
	const double axyy   = d2avector[2];
	const double ayxx   = d2avector[3];
	const double ayxy   = d2avector[4];
	const double ayyy   = d2avector[5];
	const double asxx   = d2avector[6];
	const double asxy   = d2avector[7];
	const double asyy   = d2avector[8];

	const double qdmc   = beamParams->charge / beamParams->mass / SpeedOfLight;
	const double b0g0   = beamParams->rigidity * qdmc;
	const double gamma0 = sqrt(1+b0g0*b0g0);
	const double beta0  = b0g0/gamma0;

	const double px1    = ps0[1] - ax;
	const double py1    = ps0[3] - ay;
	const double dp1    = ps0[5];

	const double h1     = sqrt( (1/beta0 + dp1)*(1/beta0 + dp1) - px1*px1 - py1*py1 - 1/b0g0/b0g0 );
	const double h3     = h1*h1*h1;

	d2H[0]           = (px1*axx + py1*ayx)*(px1*axx + py1*ayx)/h3 - (px1*axxx + py1*ayxx - axx*axx - ayx*ayx)/h1 - asxx;
	d2H[1]           = -(((1/beta0 + dp1)*(1/beta0 + dp1) - py1*py1 - 1/b0g0/b0g0)*axx + px1*py1*ayx)/h3;
	d2H[2]           = (px1*axy + py1*ayy)*(px1*axx + py1*ayx)/h3 - (px1*axxy + py1*ayxy - axx*axy - ayx*ayy)/h1 - asxy;
	d2H[3]           = -py1*(px1*axx + py1*ayx)/h3 - ayx/h1;
	d2H[4]           = 0;
	d2H[5]           = (1/beta0 + dp1)*(px1*axx + py1*ayx)/h3;

	d2H[6]           = d2H[1];
	d2H[7]           =  ( (1/beta0 + dp1)*(1/beta0 + dp1) - py1*py1 - 1/b0g0/b0g0 )/h3;
	d2H[8]           = -(((1/beta0 + dp1)*(1/beta0 + dp1) - py1*py1 - 1/b0g0/b0g0 )*axy + px1*py1*ayy)/h3;
	d2H[9]           = px1*py1/h3;
	d2H[10]          = 0;
	d2H[11]          = -(1/beta0 + dp1)*px1/h3;

	d2H[12]           = d2H[2];
	d2H[13]           = d2H[8];
	d2H[14]           = (px1*axy + py1*ayy)*(px1*axy + py1*ayy)/h3 - (px1*axyy + py1*ayyy - axy*axy - ayy*ayy)/h1 - asyy;
	d2H[15]           = -py1*(px1*axy + py1*ayy)/h3 - ayy/h1;
	d2H[16]           = 0;
	d2H[17]           = (1/beta0 + dp1)*(px1*axy + py1*ayy)/h3;

	d2H[18]           = d2H[3];
	d2H[19]           = d2H[9];
	d2H[20]           = d2H[15];
	d2H[21]           = ( (1/beta0 + dp1)*(1/beta0 + dp1) - px1*px1 - 1/b0g0/b0g0 )/h3;
	d2H[22]           = 0;
	d2H[23]           =  -(1/beta0 + dp1)*py1/h3;

	d2H[24]           = d2H[4];
	d2H[25]           = d2H[10];
	d2H[26]           = d2H[16];
	d2H[27]           = d2H[22];
	d2H[28]           = 0;
	d2H[29]           = 0;

	d2H[30]           = d2H[5];
	d2H[31]           = d2H[11];
	d2H[32]           = d2H[17];
	d2H[33]           = d2H[23];
	d2H[34]           = d2H[29];
	d2H[35]           = (px1*px1 + py1*py1 + 1/b0g0/b0g0)/h3;
}


double ffactorial(int x)
{
	double f;
	int n;

	f = (double) x;

	for(n=x-1; n>1; n--)
		f *= n;

	return f;
}


double fpower(double x, int n)
{
	double x1 = 1;
	int nn;

	for(nn=1; nn<=n; nn++)
		x1 *= x;

	return x1;
}


void Avector(double* avector, double* ctable, double* mvalues, int msize, int ssize, int lsize, int sindx, double* ps0)
{
	const double x = ps0[0];
	const double y = ps0[2];

	const double r = sqrt(x*x + y*y);
	const double theta = atan2(y,x);

	double ar = 0;

	int mindx;
	for(mindx=0; mindx<msize; mindx++)
	{
		int m = (int) mvalues[mindx];
		double cosmtheta = cos(m*theta);

		double mfact = ffactorial(m-1);
		int lfact = 1;

		int l;

		for(l=0; l<lsize/2-1; l++)
		{
			double coeff;

			if(l>0) lfact = lfact*l;

			coeff = fpower(-1,l+1) * mfact / fpower(2,2*l) / lfact / ffactorial(l+m);

			ar   += 2 * coeff * ctable[(2*l+1)*(msize*ssize) + sindx*msize + mindx] * fpower(r,2*l+m+1) * cosmtheta;
		}

	}

	avector[0] = ar * cos(theta) / beamParams->rigidity;
	avector[1] = ar * sin(theta) / beamParams->rigidity;
}


void dAvector(double* davector, double* ctable, double* mvalues, int msize, int ssize, int lsize, int sindx, double* ps0)
{
	const double x          = ps0[0];
	const double y          = ps0[2];

	const double r          = sqrt(x*x + y*y);
	const double theta      = atan2(y,x);

	const double costheta   = cos(theta);
	const double sintheta   = sin(theta);

	const double cossquared = costheta*costheta;
	const double cossin     = costheta*sintheta;
	const double sinsquared = sintheta*sintheta;

	double ardr             = 0;
	double dardr            = 0;
	double dardthetadr      = 0;
	double dasdr            = 0;
	double dasdthetadr      = 0;

	int mindx;
	for(mindx=0; mindx<msize; mindx++)
	{
		int m = (int) mvalues[mindx];

		double cosmtheta = cos(m*theta);
		double sinmtheta = sin(m*theta);

		double mfact = ffactorial(m-1);
		int lfact = 1;

		int l;

		for(l=0; l<lsize/2-1; l++)
		{
			double coeff;
			double coeffr;
			double coeffrcos;
			double coeffrsin;
			double coeffs;
			double coeffscos;
			double coeffssin;

			if(l>0) lfact = lfact*l;

			coeff        = fpower(-1,l+1) * mfact / fpower(2,2*l) / lfact / ffactorial(l+m);

			coeffr       = coeff  * ctable[(2*l+1)*(msize*ssize) + sindx*msize + mindx] * fpower(r,2*l+m);
			coeffrcos    = coeffr * cosmtheta;
			coeffrsin    = coeffr * sinmtheta;

			ardr        += 2 * coeffrcos                ;
			dardr       += 2 * coeffrcos * (2*l + m + 1);
			dardthetadr -= 2 * coeffrsin *        m     ;

			coeffs       =-coeff  * (2*l + m) * ctable[2*l*(msize*ssize) + sindx*msize + mindx] * fpower(r,2*l+m-1);
			coeffscos    = coeffs * cosmtheta;
			coeffssin    = coeffs * sinmtheta;

			dasdr       += 2 * coeffscos * (2*l + m);
			dasdthetadr -= 2 * coeffssin *        m ;
		}

	}

	davector[0] = ( ardr*sinsquared + dardr*cossquared - dardthetadr*cossin    ) / beamParams->rigidity;
	davector[1] = (-ardr*cossin     + dardr*cossin     + dardthetadr*cossquared) / beamParams->rigidity;
	davector[2] = (-ardr*cossin     + dardr*cossin     - dardthetadr*sinsquared) / beamParams->rigidity;
	davector[3] = ( ardr*cossquared + dardr*sinsquared + dardthetadr*cossin    ) / beamParams->rigidity;
	davector[4] = (                   dasdr*costheta   - dasdthetadr*sintheta  ) / beamParams->rigidity;
	davector[5] = (                   dasdr*sintheta   + dasdthetadr*costheta  ) / beamParams->rigidity;
}


void d2Avector(double* d2avector, double* ctable, double* mvalues, int msize, int ssize, int lsize, int sindx, double* ps0)
{
	const double x             = ps0[0];
	const double y             = ps0[2];

	const double r             = sqrt(x*x + y*y);
	const double theta         = atan2(y,x);

	const double costheta      = cos(theta);
	const double sintheta      = sin(theta);

	const double cossquared    = costheta*costheta;
	const double cossin        = costheta*sintheta;
	const double sinsquared    = sintheta*sintheta;

	const double coscubed      = costheta*cossquared;
	const double cossquaredsin = sintheta*cossquared;
	const double cossinsquared = costheta*sinsquared;
	const double sincubed      = sintheta*sinsquared;

	double c00;
	double c10;
	double c01;
	double c20;
	double c11;
	double c02;

	double ardr2               = 0;

	double dardrdr             = 0;
	double dardthetadr2        = 0;

	double d2ardr2             = 0;
	double d2ardthetadrdr      = 0;
	double d2ardthetadthetadr2 = 0;

	double dasdrdr             = 0;
	double dasdthetadr2        = 0;

	double d2asdr2             = 0;
	double d2asdthetadrdr      = 0;
	double d2asdthetadthetadr2 = 0;

	int mindx;

	for(mindx=0; mindx<msize; mindx++)
	{
		int m = (int) mvalues[mindx];

		double cosmtheta = cos(m*theta);
		double sinmtheta = sin(m*theta);

		double mfact = ffactorial(m-1);
		int lfact = 1;

		int l;

		for(l=0; l<lsize/2-1; l++)
		{
			double coeff;
			double coeffr;
			double coeffrcos;
			double coeffrsin;
			double coeffs;
			double coeffscos;
			double coeffssin;

			if(l>0) lfact = lfact*l;

			coeff                = fpower(-1,l+1) * mfact / fpower(2,2*l) / lfact / ffactorial(l+m);

			coeffr               = coeff  * 2 * ctable[(2*l+1)*(msize*ssize) + sindx*msize + mindx] * fpower(r,2*l+m-1);
			coeffrcos            = coeffr * cosmtheta;
			coeffrsin            = coeffr * sinmtheta;

			ardr2               += coeffrcos                            ;

			dardrdr             += coeffrcos * (2*l + m + 1)            ;
			dardthetadr2        -= coeffrsin *        m                 ;

			d2ardr2             += coeffrcos * (2*l + m + 1) * (2*l + m);
			d2ardthetadrdr      -= coeffrsin * (2*l + m + 1) *        m ;
			d2ardthetadthetadr2 -= coeffrcos *        m      *        m ;

			coeffs               =-coeff  * (2*l + m) * ctable[2*l*(msize*ssize) + sindx*msize + mindx] * fpower(r,2*l+m-2);
			coeffscos            = coeffs * cosmtheta;
			coeffssin            = coeffs * sinmtheta;

			dasdrdr             += coeffscos * (2*l + m);
			dasdthetadr2        -= coeffssin *        m ;

			d2asdr2             += coeffscos * (2*l + m) * (2*l + m - 1);
			d2asdthetadrdr      -= coeffssin * (2*l + m) *        m     ;
			d2asdthetadthetadr2 -= coeffscos *        m  *        m     ;
		}
	}

	// d^2A_x/dx^2
	c00 = -3*cossinsquared;
	c10 =  3*cossinsquared;
	c01 = -sintheta + 3*cossquaredsin - sincubed;
	c20 =  coscubed;
	c11 = -2*cossquaredsin;
	c02 =  cossinsquared;

	d2avector[0] = ( c00*ardr2 + c10*dardrdr + c01*dardthetadr2 + c20*d2ardr2 + c11*d2ardthetadrdr + c02*d2ardthetadthetadr2 ) / beamParams->rigidity;

	// d^2A_x/dx/dy
    c00 =-sintheta/4 + 9*cossquaredsin/4 - 3*sincubed/4;
    c10 = sintheta/4 - 9*cossquaredsin/4 + 3*sincubed/4;
    c01 =-coscubed   + 3*cossinsquared;
    c20 = sintheta/4 + 3*cossquaredsin/4 -   sincubed/4;
    c11 = costheta/2 +   coscubed/2      - 3*cossinsquared/2;
    c02 =-cossquaredsin;

	d2avector[1] = ( c00*ardr2 + c10*dardrdr + c01*dardthetadr2 + c20*d2ardr2 + c11*d2ardthetadrdr + c02*d2ardthetadthetadr2 ) / beamParams->rigidity;

	// d^2A_x/dy^2
   c00 =  -costheta/4    - 3*coscubed/4 + 9*cossinsquared/4;
   c10 =   costheta/4    + 3*coscubed/4 - 9*cossinsquared/4;
   c01 =-4*cossquaredsin;
   c20 =   cossinsquared;
   c11 = 2*cossquaredsin;
   c02 =   coscubed;

	d2avector[2] = ( c00*ardr2 + c10*dardrdr + c01*dardthetadr2 + c20*d2ardr2 + c11*d2ardthetadrdr + c02*d2ardthetadthetadr2 ) / beamParams->rigidity;

	// d^2A_y/dx^2
    c00 =  -sintheta/4    + 9/4*cossquaredsin - 3*sincubed/4;
    c10 =   sintheta/4    - 9/4*cossquaredsin + 3*sincubed/4;
    c01 = 4*cossinsquared;
    c20 =   cossquaredsin;
    c11 =-2*cossinsquared;
    c02 =   sincubed;

	d2avector[3] = ( c00*ardr2 + c10*dardrdr + c01*dardthetadr2 + c20*d2ardr2 + c11*d2ardthetadrdr + c02*d2ardthetadthetadr2 ) / beamParams->rigidity;

	// d^2A_y/dx/dy
    c00 =  -costheta/4    - 3*coscubed/4      + 9*cossinsquared/4;
    c10 =   costheta/4    + 3*coscubed/4      - 9*cossinsquared/4;
    c01 =-3*cossquaredsin +   sincubed;
    c20 =   costheta/4    -   coscubed/4      + 3*cossinsquared/4;
    c11 =  -sintheta/2    + 3*cossquaredsin/2 - sincubed/2;
    c02 =  -cossinsquared;

	d2avector[4] = ( c00*ardr2 + c10*dardrdr + c01*dardthetadr2 + c20*d2ardr2 + c11*d2ardthetadrdr + c02*d2ardthetadthetadr2 ) / beamParams->rigidity;

	// d^2A_y/dy^2
    c00 =-3*cossquaredsin;
    c10 = 3*cossquaredsin;
    c01 =   costheta    + coscubed   - 3*cossinsquared;
    c20 =   sincubed;
    c11 =   costheta/2  - coscubed/2 + 3*cossinsquared/2;
    c02 =   cossquaredsin;

	d2avector[5] = ( c00*ardr2 + c10*dardrdr + c01*dardthetadr2 + c20*d2ardr2 + c11*d2ardthetadrdr + c02*d2ardthetadthetadr2 ) / beamParams->rigidity;

	// d^2A_s/dx^2
    c10 =   sinsquared;
    c01 = 2*cossin;
    c20 =   cossquared;
    c11 =-2*cossin;
    c02 =   sinsquared;

	d2avector[6] = ( c10*dardrdr + c01*dardthetadr2 + c20*d2ardr2 + c11*d2ardthetadrdr + c02*d2ardthetadthetadr2 ) / beamParams->rigidity;

	// d^2A_s/dx/dy
    c10 =-cossin;
    c01 =-cossquared + sinsquared;
    c20 = cossin;
    c11 = cossquared - sinsquared;
    c02 =-cossin;

	d2avector[7] = ( c10*dardrdr + c01*dardthetadr2 + c20*d2ardr2 + c11*d2ardthetadrdr + c02*d2ardthetadthetadr2 ) / beamParams->rigidity;

	// d^2A_s/dy^2
    c10 = 3*cossquaredsin;
    c01 = costheta + coscubed - 3*cossinsquared;
    c20 = sincubed;
    c11 = costheta/2 - coscubed/2 + 3*cossinsquared/2;
    c02 = cossquaredsin;

	d2avector[8] = ( c10*dardrdr + c01*dardthetadr2 + c20*d2ardr2 + c11*d2ardthetadrdr + c02*d2ardthetadthetadr2 ) / beamParams->rigidity;
}


void TrackFieldMapGenGradients_(double* parameters)
{
	const double ax     = parameters[_apertureX];
	const double ay     = parameters[_apertureY];

	const int msize     = (int) parameters[_msize];
	const int ssize     = (int) parameters[_ssize];
	const int lsize     = (int) parameters[_lsize];

	double* mvalues     = &parameters[_mvalues];
	double* svalues     = &parameters[_svalues];
	double* ctable      = &parameters[_ctable];

	const double qdmc   = beamParams->charge / beamParams->mass / SpeedOfLight;
	const double b0g0   = beamParams->rigidity * qdmc;
	const double gamma0 = sqrt(1+b0g0*b0g0);
	const double beta0  = b0g0/gamma0;

	const double ds     = (svalues[ssize-1] - svalues[0])/(ssize-1);

	const double length =  svalues[ssize-1] - svalues[0] + ds;

	double Smatrix[6][6];

	int n;

	int row;
	int col;

	for(row=0; row<6; row++)
		for(col=0; col<6; col++)
			Smatrix[row][col] = 0;

	Smatrix[0][1] = 1;
	Smatrix[1][0] =-1;

	Smatrix[2][3] = 1;
	Smatrix[3][2] =-1;

	Smatrix[4][5] = 1;
	Smatrix[5][4] =-1;

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
			double avector[2];

			int ni;
			int ns=0;

			for(ni=0; ni<6; ni++)
				ps0[ni] = beamPS[n*6 + ni];

			Avector(avector, ctable, mvalues, msize, ssize, lsize, ns, ps0);
			ps0[1] += avector[0];
			ps0[3] += avector[1];

			for(ns=0; ns<ssize; ns++)
			{
//				double avector[2];
				double davector[6];
				double d2avector[9];

				double psm[6];

				double dH[6];
				double d2H[36];

				double k[6];
				double f[6];
				double j[36];

				int itern;

//				fprintf(logfile,"Step: %d\n",ns);

				// Calculate the vector potential and its derivatives at the initial particle co-ordinates
				Avector(avector, ctable, mvalues, msize, ssize, lsize, ns, ps0);
				dAvector(davector, ctable, mvalues, msize, ssize, lsize, ns, ps0);

				// Calculate the derivatives of the Hamiltonian at the initial phase space co-ordinates
				gradH(dH, avector, davector, ps0);

				// Make an initial estimate of the midpoint phase space co-ordinates
				for(row=0; row<6; row++)
				{
					psm[row] = ps0[row];
					for(col=0; col<6; col++)
						psm[row] += ds*Smatrix[row][col]*dH[col]/2;
				}

				// Improve the estimate... the following should be iterated until f is small enough...
				for(itern=0; itern<4; itern++)
				{
					double resid = 0;

					Avector(avector, ctable, mvalues, msize, ssize, lsize, ns, psm);
					dAvector(davector, ctable, mvalues, msize, ssize, lsize, ns, psm);
					d2Avector(d2avector, ctable, mvalues, msize, ssize, lsize, ns, psm);

					gradH(dH, avector, davector, psm);
					gradgradH(d2H, avector, davector, d2avector, psm);

					for(row=0; row<6; row++)
					{
						k[row] = 0;
						for(col=0; col<6; col++)
						{
							int nn;
							k[row] += ds*Smatrix[row][col]*dH[col]/2;
							j[6*row + col] = 0;
							if(col==row)
								j[6*row + col] = 1;
							for(nn=0; nn<6; nn++)
								j[6*row + col] -= ds*Smatrix[row][nn]*d2H[6*nn + col]/2;
						}
						f[row] = psm[row] - k[row] - ps0[row];
					}

					gaussj(j,f);

					for(row=0; row<6; row++)
					{
						psm[row] -= f[row];
						resid += f[row]*f[row];
					}

//					fprintf(logfile,"residual = %e\n",resid);

					if(resid<1e-24)
						itern = 4;
				}

				// Calculate the gradient of the Hamiltonian at the midpoint
				gradH(dH, avector, davector, psm);

				// Apply the change in phase space co-ordinates
				for(row=0; row<6; row++)
					for(col=0; col<6; col++)
						ps0[row] += ds*Smatrix[row][col]*dH[col];
			}

			Avector(avector, ctable, mvalues, msize, ssize, lsize, ssize-1, ps0);
			ps0[1] -= avector[0];
			ps0[3] -= avector[1];

			// Update the beam with the final phase space co-ordinates
			for(ni=0; ni<6; ni++)
				beamPS[n*6 + ni] = ps0[ni];

//			if(beamParams->g != 0) {
//				theta = beamSpins[n*2];
//				phi   = beamSpins[n*2 + 1];
//
//				RotateSpins(Bx0, By0, Bz0, ds, px1, py1, dp0, &theta, &phi);
//
//				beamSpins[n*2] = theta;
//				beamSpins[n*2 + 1] = phi;
//			}

			distance[n*2] += length;

		}
	}

	beamParams->globaltime += length / (beta0*SpeedOfLight);
}


void TrackFieldMapGenGradients(FieldMapGenGradientsParameters_t fieldmapgg)
{
	const int msize     = fieldmapgg.ctablemsize;
	const int ssize     = fieldmapgg.ctablessize;
	const int lsize     = fieldmapgg.ctablelsize;

	int n;

	double* parameters = (double*)calloc(5 + msize + ssize + msize*ssize*lsize,sizeof(double));

	parameters[_apertureX] = fieldmapgg.apertureX;
	parameters[_apertureY] = fieldmapgg.apertureY;

	parameters[_msize]     = (double) fieldmapgg.ctablemsize;
	parameters[_ssize]     = (double) fieldmapgg.ctablessize;
	parameters[_lsize]     = (double) fieldmapgg.ctablelsize;

	for(n=0; n<msize; n++)
		parameters[_mvalues + n] = fieldmapgg.mvalues[n];

	for(n=0; n<ssize; n++)
		parameters[_svalues + n] = fieldmapgg.svalues[n];

	for(n=0; n<msize*ssize*lsize; n++)
		parameters[_ctable + n] = fieldmapgg.ctable[n];

	TrackFieldMapGenGradients_(parameters);

	free(parameters);
}


void CopyFieldMapGenGradients(FieldMapGenGradientsParameters_t fieldmapgg)
{
	const int msize     = fieldmapgg.ctablemsize;
	const int ssize     = fieldmapgg.ctablessize;
	const int lsize     = fieldmapgg.ctablelsize;

	int n;

	double* parameters = (double*)calloc(5 + msize + ssize + msize*ssize*lsize,sizeof(double));

	parameters[_apertureX] = fieldmapgg.apertureX;
	parameters[_apertureY] = fieldmapgg.apertureY;

	parameters[_msize]     = (double) fieldmapgg.ctablemsize;
	parameters[_ssize]     = (double) fieldmapgg.ctablessize;
	parameters[_lsize]     = (double) fieldmapgg.ctablelsize;

	for(n=0; n<msize; n++)
		parameters[_mvalues + n] = fieldmapgg.mvalues[n];

	for(n=0; n<ssize; n++)
		parameters[_svalues + n] = fieldmapgg.svalues[n];

	for(n=0; n<msize*ssize*lsize; n++)
		parameters[_ctable + n] = fieldmapgg.ctable[n];

	AppendComponent(TrackFieldMapGenGradients_, parameters);
}

