#include <math.h>
#include <stdlib.h>

#include "CParticleTracking.h"
#include "GlobalVariables.h"
#include "AppendComponent.h"
//#include "TrackingLog.h"

void ParticleTracking(InitialisationParameters_t* params)
{
	params->valid = 1;
	params->version = "Particle tracking in C version 3.1 (11 August 2015)";
}

double Precision()
{
	return (double)1e-8;
}

void InitialiseTracking(BeamParameters_t* beamParams_, double* beamPS_, double* beamSpins_, double* distance_)
{
	beamParams = beamParams_;
	beamPS = beamPS_;
	beamSpins = beamSpins_;
	distance = distance_;

//	logfile = fopen("D:/Documents/Research/Projects/HL-LHC/Fringe Fields/RK Integrator/logfile.dat","w");
}

void TrackBeamline(int n1, int n2)
{
	BeamlineComponent_t* blcpt = beamline;

	int n = 1;

	while(blcpt && (n<n1))
	{
		blcpt = blcpt->nextComponent;
		n++;
	}

	while(blcpt && (n<=n2))
	{
		(blcpt->trackingRoutine)(blcpt->parameterList);
		blcpt = blcpt->nextComponent;
		n++;
	}
}

void DeleteBeamlineCopy()
{
	BeamlineComponent_t* blcpt1 = beamline;
	BeamlineComponent_t* blcpt2 = beamline;

	while(blcpt1)
	{
		free(blcpt1->parameterList);
		blcpt2 = blcpt1->nextComponent;
		free(blcpt1);
		blcpt1 = blcpt2;
	}

	beamline = NULL;
}

void FinishTracking()
{
//	fclose(logfile);
}

void RotateSpins(double bx, double by, double bz, double ds, double px, double py, double dp, double* theta, double* phi)
{
	double polsnx0;
	double polsny0;
	double polsnz0;

	double bgRef;
	double gammaRef;
	double betaRef;

	double gamma;
	double beta;

	double ptot;
	double pz2;
	double pz;

	double pdotb;

	double bParx;
	double bPary;
	double bParz;

	double bPerpx;
	double bPerpy;
	double bPerpz;

	double emdg;
	double bigG;

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

	polsnx0  = sin(*theta)*cos(*phi);
	polsny0  = sin(*theta)*sin(*phi);
	polsnz0  = cos(*theta);

	bgRef    = beamParams->rigidity * beamParams->charge / beamParams->mass / SpeedOfLight;
	gammaRef = sqrt(1+bgRef*bgRef);
	betaRef  = bgRef/gammaRef;

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

	emdg     = beamParams->charge / beamParams->mass / gamma;
	bigG     = (beamParams->g - 2)/2;
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

	*theta = acos(polsnz1);
	*phi   = atan2(polsny1,polsnx1);
}
