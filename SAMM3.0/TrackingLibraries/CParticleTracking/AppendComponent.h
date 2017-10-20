typedef struct CParticleTracking_BeamlineComponent
{
	void (*trackingRoutine)(double*);
	double* parameterList;
	struct CParticleTracking_BeamlineComponent* nextComponent;
} BeamlineComponent_t;

void AppendComponent(void (*trackingRoutine)(double*), double* parameters);
