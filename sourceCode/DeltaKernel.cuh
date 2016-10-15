#include "DeltaPreviousPoint.cuh"


#include<strstream>

__device__ double dfabs(double b1, double b2);
__device__ int ifabs(int b1, int b2);
__device__ int findAminoInPrior(char *amino, char aa, int size);
__device__ int findPtmIndex(int* numptm, int t);
__host__ void convert(vector<mass_t>& mass, vector<int>& massInt, vector<mass_t>& peaks, vector<int>& peaksInt, vector<mass_t>& delta, vector<int>& deltaInt);
__host__ void calcDiff(vector<int>& massInt, vector<int>& peaksInt, vector<int>& diff);
__host__ void calcDiff(vector<double>& mass, vector<double>& peaks, vector<double>& diff);
__global__ void KernelInitialDeltaDiag2(int x, int y, int* diff,int* numDeltaCoords, int* deltaCoordsx, int* deltaCoordsy, int tolerance, int* delta, int deltaSize);
__global__ void KernelInitialDeltaDiag(int x, int y, double* diff, int* numDeltaCoords, int* deltaCoordsx, int* deltaCoordsy, double tolerance, double* delta, int deltaSize);