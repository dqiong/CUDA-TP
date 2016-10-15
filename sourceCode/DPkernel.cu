#include "DPkernel.cuh"

__host__ __device__ int maxInt(int a, int b){
	if (a > b)
		return a;
	else return b;
}
__host__ void cpu_dp_kernel(int x, int y, thrust::host_vector<int>& m0, thrust::host_vector<int>& d0, thrust::host_vector<int>& m1, thrust::host_vector<int>& d1, thrust::host_vector<int>& previousMatchx, thrust::host_vector<int>& previousMatchy, thrust::host_vector<int>& numDeltaDiag, thrust::host_vector<int>& deltaCoordsx, thrust::host_vector<int>& deltaCoordsy){

	//clock_t start=clock();

	for (int j = 0; j < y; j++){
		for (int i = 0; i < x; i++){

			int index = j*x + i;
			if (index != 0){
				int diag = d1[previousMatchy[index] * x + previousMatchx[index]];
				//if (diag == 0)
				//	diag = -1;
				int mvalue = 0;
				if (i>0 && j>0){
					mvalue = m0[(j - 1)*x + i - 1];
				}
				int diagDelta = 0;

				if (numDeltaDiag[index] > 0){
					int max = 0;
					for (int k = 0; k < 5; k++){
						int deltaIndexi = deltaCoordsx[5 * index + k];
						int deltaIndexj = deltaCoordsy[5 * index + k];
						int deltaIndex = deltaIndexj*x + deltaIndexi;
						int deltatemp = d0[deltaIndex];
						if (deltatemp>max)
							max = deltatemp;
					}
					diagDelta = max;
				}
				d1[index] = maxInt(maxInt(diag + 1, mvalue), diagDelta + 1);

				//std::cout << d1[index] << endl;

				int m_left = 0;
				int m_up = 0;
				if (i > 0)
					m_left = m1[j*x + i - 1];
				if (j > 0)
					m_up = m1[(j - 1)*x + i];
				m1[index] = maxInt(maxInt(d1[index], m_left), m_up);
			}

		}


	}
	//clock_t end=clock();
	//std::cout << "CPU 计算一次M和D表:" << end - start << "ms" << endl;
}


__global__ void gpu_dp_kernel(int x, int y, int* m0, int* d0, int* m1, int* d1, int* previousMatchx, int* previousMatchy, int* numDeltaDiag, int* deltaCoordsx, int* deltaCoordsy){

	int idx = blockDim.x*blockIdx.x + threadIdx.x;
	for (int k = 0; k < x + y - 1; k++){

		int i = k - idx;
		int j = idx;
		int index = j*x + i;
		if (index != 0 && j<y && (i>-1) && (i<x)){
			int diag = d1[previousMatchy[index] * x + previousMatchx[index]];
			//if (diag == 0)
			//	diag = -1;
			int mvalue = 0;
			if (i>0 && j>0){
				mvalue = m0[(j - 1)*x + i - 1];
			}
			int diagDelta = 0;

			if (numDeltaDiag[index] > 0){
				int max = 0;
				for (int k = 0; k < 5; k++){
					int deltaIndexi = deltaCoordsx[5 * index + k];
					int deltaIndexj = deltaCoordsy[5 * index + k];
					int deltaIndex = deltaIndexj*x + deltaIndexi;
					int deltatemp = d0[deltaIndex];
					if (deltatemp>max)
						max = deltatemp;
				}
				diagDelta = max;
			}
			d1[index] = maxInt(maxInt(diag + 1, mvalue), diagDelta + 1);

			//__syncthreads();
			//std::cout << d1[index] << endl;

			int m_left = 0;
			int m_up = 0;
			if (i > 0)
				m_left = m1[j*x + i - 1];
			if (j > 0)
				m_up = m1[(j - 1)*x + i];
			m1[index] = maxInt(maxInt(d1[index], m_left), m_up);
		}

		__syncthreads();


	}

}

__global__ void gpu_dp_kernel2(int x, int y, int* m0, int* d0, int* m1, int* d1, int* previousMatchx, int* previousMatchy, int* numDeltaDiag, int* deltaCoordsx, int* deltaCoordsy, int internalRow, int internalColumn, int rows, int remainder, int columns){
	int idx = blockDim.x*blockIdx.x + threadIdx.x;
	for (int k = 0; k < rows + columns + 3; k++){
		for (int interi = 0; interi < internalColumn; interi++){
			int i = k *internalRow - idx + interi;
			int j = idx;
			int index = j*x + i;
			if (index != 0 && i<x && (i>-1)){
				/*
				if (i<0){
				i = i - remainder;
				j = j - internalColumn + 1;
				}
				*/
				int diag = d1[previousMatchy[index] * x + previousMatchx[index]];
				//if (diag == 0)
				//	diag = -1;
				int mvalue = 0;
				if (i > 0 && j > 0){
					mvalue = m0[(j - 1)*x + i - 1];
				}
				int diagDelta = 0;

				if (numDeltaDiag[index] > 0){
					int max = 0;
					for (int k = 0; k < 5; k++){
						int deltaIndexi = deltaCoordsx[5 * index + k];
						int deltaIndexj = deltaCoordsy[5 * index + k];
						int deltaIndex = deltaIndexj*x + deltaIndexi;
						int deltatemp = d0[deltaIndex];
						if (deltatemp>max)
							max = deltatemp;
					}
					diagDelta = max;
				}
				d1[index] = maxInt(maxInt(diag + 1, mvalue), diagDelta + 1);
				//std::cout << d1[index] << endl;

				int m_left = 0;
				int m_up = 0;
				if (i > 0)
					m_left = m1[j*x + i - 1];
				if (j > 0)
					m_up = m1[(j - 1)*x + i];
				m1[index] = maxInt(maxInt(d1[index], m_left), m_up);
				

			}
			__syncthreads();
		}

	}
}
__global__ void gpu_dp_residue(int x, int y, int* m0, int* d0, int* m1, int* d1, int* previousMatchx, int* previousMatchy, int* numDeltaDiag, int* deltaCoordsx, int* deltaCoordsy, int internalRow, int internalColumn, int rows, int columns, int remainderx, int remaindery){

}

void cpu_residue(int x, int y, thrust::device_vector<int>& m0, thrust::device_vector<int>& d0, thrust::device_vector<int>& m1, thrust::device_vector<int>& d1, thrust::device_vector<int>& previousMatchx, thrust::device_vector<int>& previousMatchy, thrust::device_vector<int>& numDeltaDiag, thrust::device_vector<int>& deltaCoordsx, thrust::device_vector<int>& deltaCoordsy, int internalRow, int internalColumn, int rows, int columns, int remainderx, int remaindery){

	clock_t start = clock();
	if (remainderx != 0){
		for (int j = 0; j < y; j++){
			for (int i = x - remainderx; i < x; i++){

				int index = j*x + i;
				if (index != 0){
					int diag = d1[previousMatchy[index] * x + previousMatchx[index]];
					//if (diag == 0)
					//	diag = -1;
					int mvalue = 0;
					if (i>0 && j>0){
						mvalue = m0[(j - 1)*x + i - 1];
					}
					int diagDelta = 0;

					if (numDeltaDiag[index] > 0){
						int max = 0;
						for (int k = 0; k < 5; k++){
							int deltaIndexi = deltaCoordsx[5 * index + k];
							int deltaIndexj = deltaCoordsy[5 * index + k];
							int deltaIndex = deltaIndexj*x + deltaIndexi;
							int deltatemp = d0[deltaIndex];
							if (deltatemp>max)
								max = deltatemp;
						}
						diagDelta = max;
					}
					d1[index] = maxInt(maxInt(diag + 1, mvalue), diagDelta + 1);

					//std::cout << d1[index] << endl;

					int m_left = 0;
					int m_up = 0;
					if (i > 0)
						m_left = m1[j*x + i - 1];
					if (j > 0)
						m_up = m1[(j - 1)*x + i];
					m1[index] = maxInt(maxInt(d1[index], m_left), m_up);
				}

			}


		}
	}
	if (remaindery != 0){
		for (int j = y - remaindery; j < y; j++){
			for (int i = 0; i < x; i++){

				int index = j*x + i;
				if (index != 0){
					int diag = d1[previousMatchy[index] * x + previousMatchx[index]];
					//if (diag == 0)
					//	diag = -1;
					int mvalue = 0;
					if (i>0 && j>0){
						mvalue = m0[(j - 1)*x + i - 1];
					}
					int diagDelta = 0;

					if (numDeltaDiag[index] > 0){
						int max = 0;
						for (int k = 0; k < 5; k++){
							int deltaIndexi = deltaCoordsx[5 * index + k];
							int deltaIndexj = deltaCoordsy[5 * index + k];
							int deltaIndex = deltaIndexj*x + deltaIndexi;
							int deltatemp = d0[deltaIndex];
							if (deltatemp>max)
								max = deltatemp;
						}
						diagDelta = max;
					}
					d1[index] = maxInt(maxInt(diag + 1, mvalue), diagDelta + 1);

					//std::cout << d1[index] << endl;

					int m_left = 0;
					int m_up = 0;
					if (i > 0)
						m_left = m1[j*x + i - 1];
					if (j > 0)
						m_up = m1[(j - 1)*x + i];
					m1[index] = maxInt(maxInt(d1[index], m_left), m_up);
				}

			}


		}
	}
	clock_t end = clock();
	cout << "residue:" << end - start << endl;

}