
#include <math.h>

#include "CParticleTracking.h"

BeamParameters_t* beamParams;
double* beamPS;
double* distance;

void TrackBPM(BPMParameters_t bpm, double* bpmX, double* bpmY, int* nparticles)
{
	const double ax = bpm.apertureX;
	const double ay = bpm.apertureY;

	int n;

	*bpmX = 0;
	*bpmY = 0;
	*nparticles = 0;

	for(n=0; n<beamParams->nparticles; n++)
	{
		if(distance[n*2 + 1])
		{
			double x0  = beamPS[n*6];
			double y0  = beamPS[n*6 + 2];

			(*nparticles)++;
			*bpmX += x0;
			*bpmY += y0;

			if( (ax>0) & (ay>0) )
				if( (x0/ax)*(x0/ax) + (y0/ay)*(y0/ay) > 1)
					distance[n*2 + 1] = 0;
		}
	}

	if(*nparticles>0)
	{
		*bpmX /= *nparticles;
		*bpmY /= *nparticles;
	}
}
