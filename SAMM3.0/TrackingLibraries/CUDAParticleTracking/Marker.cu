#include <cuda_runtime.h>

#include "CUDAParticleTracking.h"
#include "TrackCUDA.h"

__global__ void MarkerMap(Component* drift, Beam* beam, Particle* ptcle)
{
	// A marker has no effect on the beam
}

extern "C" __host__ void TrackMarkerCUDA()
{
	// A marker has no effect on the beam
}

extern "C" __host__ void CopyMarkerCUDA()
{
	AppendComponent(&MarkerMap, NULL);
}
