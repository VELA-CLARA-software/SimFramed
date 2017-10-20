
#define threadsPerBlock 128

#define _nbeamparams 7
#define _globaltime 0
#define _mass 1
#define _charge 2
#define _g 3
#define _rigidity 4
#define _gamma 5
#define _beta 6

#define _nptcleparams 10
#define _x(index)     (_nptcleparams*(index))
#define _px(index)    (_nptcleparams*(index)+1)
#define _y(index)     (_nptcleparams*(index)+2)
#define _py(index)    (_nptcleparams*(index)+3)
#define _ct(index)    (_nptcleparams*(index)+4)
#define _dp(index)    (_nptcleparams*(index)+5)
#define _theta(index) (_nptcleparams*(index)+6)
#define _phi(index)   (_nptcleparams*(index)+7)
#define _s(index)     (_nptcleparams*(index)+8)
#define _f(index)     (_nptcleparams*(index)+9)

typedef double Particle;

typedef double Component;

typedef double Beam;

typedef struct CUDAParticleTracking_BeamlineComponent
{
	void (*trackingRoutine)(Component*, Beam*, Particle*);
	double* parameterList;
	struct CUDAParticleTracking_BeamlineComponent* nextComponent;
} BeamlineComponent_t;

extern "C" __host__ void SetBeam();

extern "C" __host__ void GetBeam();

extern "C" __host__ void AppendComponent(void (*trackingRoutine)(Component*, Beam*, Particle*), Component* parameterList);
