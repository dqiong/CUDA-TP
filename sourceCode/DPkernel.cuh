#include <cuda.h>
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <cuda.h>
#include <device_functions.h>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>

#include <thrust/copy.h>
#include <thrust/sort.h>

#include<time.h>
#include<iostream>
#include<vector>
using namespace std;


__host__ void cpu_dp_kernel(int x, int y, thrust::host_vector<int>& m0, thrust::host_vector<int>& d0, thrust::host_vector<int>& m1, thrust::host_vector<int>& d1, thrust::host_vector<int>& previousMatchx, thrust::host_vector<int>& previousMatchy, thrust::host_vector<int>& numDeltaDiag, thrust::host_vector<int>& deltaCoordsx, thrust::host_vector<int>& deltaCoordsy);

__global__ void gpu_dp_kernel(int x,int y,int* m0,int* d0,int* m1,int* d1,int* previousMatchx,int* previousMatchy,int* numDeltaDiag,int* deltaCoordsx,int* deltaCoordsy);

__global__ void gpu_dp_kernel2(int x,int y,int* m0,int* d0,int* m1,int* d1,int* previousMatchx,int* previousMatchy,int* numDeltaDiag,int* deltaCoordsx,int* deltaCoordsy,int internalRow,int internalColumn,int rows,int remainder,int columns);

__global__ void gpu_dp_residue(int x, int y, int* m0, int* d0, int* m1, int* d1, int* previousMatchx, int* previousMatchy, int* numDeltaDiag, int* deltaCoordsx, int* deltaCoordsy, int internalRow, int internalColumn, int rows, int columns,int remainderx,int remaindery);

void cpu_residue(int x, int y,thrust::device_vector<int>& m0, thrust::device_vector<int>& d0, thrust::device_vector<int>& m1, thrust::device_vector<int>& d1, thrust::device_vector<int>& previousMatchx, thrust::device_vector<int>& previousMatchy, thrust::device_vector<int>& numDeltaDiag, thrust::device_vector<int>& deltaCoordsx, thrust::device_vector<int>& deltaCoordsy,int internalRow,int internalColumn,int rows,int columns,int remainderx,int remaindery);