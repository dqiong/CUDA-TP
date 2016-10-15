#include "DeltaKernel.cuh"

__device__ double dfabs(double b1, double b2){
	double result = b1 - b2;
	if (result > 0)
		return result;
	else return -result;
}
__device__ int ifabs(int b1, int b2){
	int result = b1 - b2;
	if (result > 0)
		return result;
	else return -result;
}
__device__ int findAminoInPrior(char *amino, char aa, int size){
	for (int i = 0; i < size; i++){
		if (amino[i] == aa)
			return i;
	}
	return -1;
}
__device__ int findPtmIndex(int* numptm, int t){
	int num = 0;
	for (int i = 0; i < t; i++){
		num += numptm[i];
	}
	return num;
}
__host__ void convert(vector<mass_t>& mass, vector<int>& massInt, vector<mass_t>& peaks, vector<int>& peaksInt, vector<mass_t>& delta, vector<int>& deltaInt){
	int fold = 10000;
	for (int i = 0; i < mass.size(); i++){
		massInt.push_back(mass[i] * fold);
	}
	for (int i = 0; i < peaks.size(); i++){
		peaksInt.push_back(peaks[i] * fold);
	}
	for (int i = 0; i < delta.size(); i++){
		deltaInt.push_back(delta[i] * fold);
	}
}
__host__ void calcDiff(vector<int>& massInt, vector<int>& peaksInt, vector<int>& diff){
	int x = massInt.size();
	int y = peaksInt.size();
	for (int j = 0; j < y; j++){
		for (int i = 0; i < x; i++){
			int index = j*x + i;
			diff[index] = peaksInt[j] - massInt[i];
		}
	}
}
__host__ void calcDiff(vector<double>& mass, vector<double>& peaks, vector<double>& diff){
	int x = mass.size();
	int y = peaks.size();
	for (int j = 0; j < y; j++){
		for (int i = 0; i < x; i++){
			int index = j*x + i;
			diff[index] = peaks[j] - mass[i];
		}
	}
}
__global__ void KernelInitialDeltaDiag2(int x, int y, int* diff, int* numDeltaCoords, int* deltaCoordsx, int* deltaCoordsy, int tolerance, int* delta, int deltaSize){
	int idx = blockDim.x*blockIdx.x + threadIdx.x;
	//atomicAdd(&num, 1);
	int numDelta = 0;
	if (idx < x*y){
		int currx = idx%x;
		int curry = idx / x;
		for (int j = 0; j < curry; j++){
			for (int i = 0; i < currx; i++){
				int index = j*x + i;
				int valueDelta = diff[idx] - diff[index];
				//int ptmIndex = findPtmIndex(numptm, faa);
				if (valueDelta>0 && valueDelta < 800000){
					for (int no_ptmi = 0; no_ptmi < deltaSize; no_ptmi++){
						int toleranceDelta = ifabs(valueDelta, delta[no_ptmi]);
						if (toleranceDelta < tolerance){
							deltaCoordsx[5 * idx + numDelta % 5] = i;
							deltaCoordsy[5 * idx + numDelta % 5] = j;
							numDelta++;
							break;
						}
					}
				}

			}
			//}

		}
		numDeltaCoords[idx] = numDelta;
	}
}
__global__ void KernelInitialDeltaDiag(int x, int y, double* diff, int* numDeltaCoords, int* deltaCoordsx, int* deltaCoordsy, double tolerance, double* delta, int deltaSize){
	int idx = blockDim.x*blockIdx.x + threadIdx.x;
	//atomicAdd(&num, 1);
	int numDelta = 0;
	if (idx < x*y){
		int currx = idx%x;
		int curry = idx / x;
		for (int j = 0; j < curry; j++){
			for (int i = 0; i < currx; i++){
				int index = j*x + i;
				int valueDelta = diff[idx] - diff[index];
				//int ptmIndex = findPtmIndex(numptm, faa);
				if (valueDelta>0 && valueDelta < 80){
					for (int no_ptmi = 0; no_ptmi < deltaSize; no_ptmi++){
						int toleranceDelta = dfabs(valueDelta, delta[no_ptmi]);
						if (toleranceDelta < tolerance){
							deltaCoordsx[5 * idx + numDelta % 5] = i;
							deltaCoordsy[5 * idx + numDelta % 5] = j;
							numDelta++;
							break;
						}
					}
				}

			}
			numDeltaCoords[idx] = numDelta;
		}
	}
}