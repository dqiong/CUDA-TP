

#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <cuda.h>
#include <device_functions.h>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/sort.h>
#include <thrust/generate.h>
#include <thrust/copy.h>

#include "PreviousPoint.h"
#include "ProteinSequence.h"
#include "MonoSpectrum.h"
#include "MonoFragments.h"

#include <iostream>
#include <math.h>
#include <vector>
#include <map>
#include <string>
#include <time.h>
#include <stdio.h>
using namespace std;


 void InitialDeltaPreviousPoint(string protein,Config config,thrust::host_vector<mass_t>& cutMasses,thrust::host_vector<mass_t>& cutMonoPeaks,int x,int y,double tolerance,thrust::host_vector<int>& hNumDeltaCoords,thrust::host_vector<int>& hDeltaCoordsx,thrust::host_vector<int>& hDeltaCoordsy);

