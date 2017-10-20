#include <stdlib.h>

#include "AppendComponent.h"
#include "GlobalVariables.h"

void AppendComponent(void (*trackingRoutine)(double*), double* parameters)
{
	static BeamlineComponent_t* currentBeamlineComponent;

	BeamlineComponent_t* blcpt = (BeamlineComponent_t *)malloc(sizeof(BeamlineComponent_t));

	blcpt->trackingRoutine = trackingRoutine;
	blcpt->parameterList   = parameters;
	blcpt->nextComponent   = NULL;

	if(beamline)
	{
		currentBeamlineComponent->nextComponent = blcpt;
		currentBeamlineComponent = blcpt;
	}
	else
	{
		beamline = blcpt;
		currentBeamlineComponent = blcpt;
	}
}
