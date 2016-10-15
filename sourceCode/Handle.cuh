
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <cuda.h>
#include <device_functions.h>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/sort.h>
#include <thrust/generate.h>
#include <thrust/copy.h>

#include <stdio.h>
#include "includes.h"

#include "MonoSpectrum.h"
#include "MonoFragments.h"
#include "ProteinSequence.h"
#include "PreviousPoint.h"
#include "DeltaPreviousPoint.cuh"
#include "DPkernel.cuh"
#include "DeltaKernel.cuh"
#include "FilterProtein.cuh"


#include<map>
#include<strstream>
#include<string>
#include<vector>
#include<iomanip>

using namespace std;

//when there is no mass shift,the proceeds of array D and M
void InitialArray(thrust::host_vector<int>& d, thrust::host_vector<int>& m, int x, int y, thrust::host_vector<CutPeak> peaks, thrust::host_vector<mass_t> masses, double tolerance);

void InitialArray(vector<int>& d, vector<int>& m, int x, int y, vector<CutPeak> peaks, vector<mass_t> masses, double tolerance);

//compute the previous point for all points in the grid
void InitialPreviousPoint(Node* &root, thrust::host_vector<mass_t>& A, thrust::host_vector<mass_t>& B, int x, int y, Point *match);
void InitialPreviousPoint(Node* &root, vector<mass_t>& A, vector<mass_t>& B, int x, int y, Point *match);

void InitialMapAminoPtm(map<char, bool>& m);

void SplitString(const string& s, vector<std::string>& v, const std::string& c);

void gpuFunction(vector<vector<mass_t>>& cutMasses,vector<CutPeak>& cutPeaks,vector<mass_t>& cutPeaksMonoMass,vector<mass_t> ptm,vector<int>& tscore);
void gpuFunction2();
void approxConvolutionmy(const vector<double>& masses, const vector<double>& peaks, double e);
__global__ void testt(int m);
__device__ bool findVal(int* a, int value,int m);
__device__ void BubbleSort(double *arr, int n);
__global__ void gpuFilter(int* cutMass_begin, int* cutMass_end, double* a,double* b,int m,int proteins);
__global__ void gpuApproxConvolution(int* cutMass_begin, int* cutMass_end, double* a, double* b, int m, int proteins);
__global__ void gpuExactConvolution(int* cutMass_begin, int* cutMass_end, double* a, double* b, int m, int proteins);
__global__ void gpuFilter();
__global__ void calScore(double* a, double* b, int n, int m, double ma, double mb, double e, int* c, int num, int* score,int max);
void filterProtein(vector<vector<mass_t>>& cutMasses, vector<MonoSpectrum>& spectrum, vector<vector<mass_t>>& cutPeaksMonoMass,vector<vector<CutPeak>>& cutPeaks);