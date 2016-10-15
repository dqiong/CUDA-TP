/*
#ifdef __CUDACC__
#define KERNEL_ARGS2(grid, block) <<< grid, block >>>
#define KERNEL_ARGS3(grid, block, sh_mem) <<< grid, block, sh_mem >>>
#define KERNEL_ARGS4(grid, block, sh_mem, stream) <<< grid, block, sh_mem, stream >>>
#else
#define KERNEL_ARGS2(grid, block)
#define KERNEL_ARGS3(grid, block, sh_mem)
#define KERNEL_ARGS4(grid, block, sh_mem, stream)
#endif
*/

#include "DeltaPreviousPoint.cuh"
#include "Handle.cuh"

#include<strstream>



void cpuKernel2(int x, int y, vector<int>& diff, vector<int>& coords, vector<int>& numDeltaCoords, vector<int>& deltaCoordsx, vector<int>& deltaCoordsy, int tolerance, vector<int>& delta, int deltaSize){
	clock_t start = clock();
	for (int cury = 0; cury < y; cury++){
		for (int curx = 0; curx < x; curx++){
			int idx = cury*x + curx;
			int numDelta = 0;
			for (int j = 0; j < cury; j++){
				for (int i = 0; i < curx; i++){
					int index = j*x + i;
					int valueDelta = diff[idx] - diff[index];
					//int ptmIndex = findPtmIndex(numptm, faa);
					if (valueDelta>0 && valueDelta < 800000){
						for (int no_ptmi = 0; no_ptmi < deltaSize; no_ptmi++){
							int toleranceDelta = abs((int)valueDelta - delta[no_ptmi]);
							if (toleranceDelta < tolerance){
								deltaCoordsx[5 * idx + numDelta % 5] = i;
								deltaCoordsy[5 * idx + numDelta % 5] = j;
								numDelta++;
								break;
							}
						}
					}
				}
			}
			numDeltaCoords[idx] = numDelta;
		}
	}
	clock_t end = clock();
	cout << "diag_delta串行时间:" << end - start << "ms" << endl;

}

/*
__global__ void t_KernelInitialDeltaDiag2(int x, int y, int* diff, int* coords, int* numDeltaCoords, int* deltaCoordsx, int* deltaCoordsy, int tolerance, int* delta, int deltaSize){
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
__global__ void t_KernelInitialDeltaDiagInt(int x, int y, int* mass, int* peaks, int* coords, int* numDeltaCoords, int* deltaCoordsx, int* deltaCoordsy, int tolerance, int* delta, int deltaSize){
	int idx = blockDim.x*blockIdx.x + threadIdx.x;
	//atomicAdd(&num, 1);
	int numDelta = 0;
	if (idx < x*y){
		int currx = idx%x;
		int curry = idx / x;
		for (int j = 0; j < curry; j++){
			for (int i = 0; i < currx; i++){
				int index = j*x + i;
				int valueCurr = peaks[curry] - mass[currx];
				int valuePrevious = peaks[j] - mass[i];
				int valueDelta = valueCurr - valuePrevious;
				//int ptmIndex = findPtmIndex(numptm, faa);
				for (int no_ptmi = 0; no_ptmi < deltaSize; no_ptmi++){
					int toleranceDelta = ifabs(valueDelta, delta[no_ptmi]);
					if (toleranceDelta < tolerance){
						deltaCoordsx[5 * idx + numDelta % 5] = i;
						deltaCoordsy[5 * idx + numDelta % 5] = j;
						numDelta++;
					}
				}
			}
			//}

		}
		numDeltaCoords[idx] = numDelta;
	}

	//__syncthreads();
}
__global__ void t_KernelInitialDeltaDiag(int x, int y, double* mass, double* peaks, int* coords, int* numDeltaCoords, int* deltaCoordsx, int* deltaCoordsy, double tolerance, double* delta, int deltaSize){
	int idx = blockDim.x*blockIdx.x + threadIdx.x;
	//atomicAdd(&num, 1);
	int numDelta = 0;
	if (idx < x*y){
		int currx = idx%x;
		int curry = idx / x;
		for (int j = 0; j < curry; j++){
			for (int i = 0; i < currx; i++){
				int index = j*x + i;
				double valueCurr = mass[currx] - peaks[curry];
				double valuePrevious = mass[i] - peaks[j];
				double valueDelta = dfabs(valueCurr, valuePrevious);
				//int ptmIndex = findPtmIndex(numptm, faa);
				for (int no_ptmi = 0; no_ptmi < deltaSize; no_ptmi++){
					double toleranceDelta = dfabs(valueDelta, delta[no_ptmi]);
					if (toleranceDelta < tolerance){
						deltaCoordsx[5 * idx + numDelta % 5] = i;
						deltaCoordsy[5 * idx + numDelta % 5] = j;
						numDelta++;
					}
				}
			}
			//}

		}
		numDeltaCoords[idx] = numDelta;
	}

	//__syncthreads();
}
*/
void InitialDeltaPreviousPoint(string protein, Config config, thrust::host_vector<mass_t>& cutMasses, thrust::host_vector<mass_t>& cutMonoPeaks, int x, int y, double tolerance, thrust::host_vector<int>& hNumDeltaCoords, thrust::host_vector<int>& hDeltaCoordsx, thrust::host_vector<int>& hDeltaCoordsy){

	//char *diagPath = "output\\InitialDiag2.txt";
	//ofstream out(diagPath);
	//out.setf(ios::showpoint);
	//out.precision(6);
	//out.setf(ios::fixed);

	//deal with amino and its PTMs
	/*
	map<char, vector<mass_t>> char2ptmMass;
	map<char, bool> aminoPtm;
	vector<mass_t> vdelta;
	InitialMapAminoPtm(aminoPtm);
	vector<PTM> list = config.getAllPTMs();
	for (int j = 0; j < list.size(); j++){
		vdelta.push_back(list[j].delta);
		//out << "ptm:" << list[j].label <<" mass:"<<list[j].delta<< " aa:";
		for (int k = 0; k < list[j].applicableAAs.size(); k++){
			char aa = config.get_aa2label()[list[j].applicableAAs[k]][0];
			vector<mass_t> ptmValue;
			if (char2ptmMass.count(aa) == 0){
				aminoPtm[aa] = true;
				for (int tempi = 0; tempi < list.size(); tempi++){
					int flag = 0;
					for (int tempk = 0; tempk < list[tempi].applicableAAs.size(); tempk++){
						if (aa == config.get_aa2label()[list[tempi].applicableAAs[tempk]][0])
							flag = 1;
					}
					if (flag == 1)
						ptmValue.push_back(list[tempi].delta);
				}
				char2ptmMass.insert(make_pair(aa, ptmValue));

			}

		}
		//out << config.get_aa2label()[list[j].applicableAAs[k]] << " ";
		//out << endl;
	}
	const char* tempSequence = protein.c_str();
	vector<char> vsequence(protein.length());
	for (int i = 0; i < protein.length(); i++){
		vsequence.push_back(tempSequence[i]);
	}
	thrust::device_vector<char> dprotein = vsequence; //protein sequence
	thrust::device_vector<mass_t> dPTM;           //the prior ptm
	thrust::device_vector<int> dNumPTM;           //the number ptm of a amino
	thrust::host_vector<int> hNumPTM;
	thrust::host_vector<mass_t> hPTM;
	thrust::device_vector<char> dAmino;           //the prior knowledge amino
	thrust::host_vector<char> hAmino;
	thrust::device_vector<double> dDelta = vdelta;
	int deltaSize = vdelta.size();

	thrust::device_vector<double> dmass = cutMasses;
	thrust::device_vector<double> dpeaks = cutMonoPeaks;
	thrust::device_vector<int> dCoords(2 * x*y);  //record result
	thrust::device_vector<int> dNumDeltaCoords(x*y);
	//thrust::host_vector<int> hNumDeltaCoords(x*y);
	thrust::device_vector<int> dDeltaCoordsx(5 * x*y, 0);
	//thrust::host_vector<int> hDeltaCoordsx(5*x*y,0);
	thrust::device_vector<int> dDeltaCoordsy(5 * x*y, 0);
	//thrust::host_vector<int> hDeltaCoordsy(5*x*y,0);

	vector<int> coords(2 * x*y);
	vector<int> numDeltaCoords(x*y);
	vector<int> coordsx(5 * x*y, 0);
	vector<int> coordsy(5 * x*y, 0);
	vector<int> massInt;
	vector<int> peaksInt;
	vector<int> deltaInt;
	vector<int> diff(x*y);
	convert(cutMasses, massInt, cutMonoPeaks, peaksInt, vdelta, deltaInt);
	calcDiff(massInt, peaksInt, diff);
	thrust::device_vector<int> dmassInt = massInt;
	thrust::device_vector<int> dpeaksInt = peaksInt;
	thrust::device_vector<int> ddeltaInt = deltaInt;
	thrust::device_vector<int> ddiff = diff;

	map<char, vector<mass_t>>::iterator it = char2ptmMass.begin();
	for (; it != char2ptmMass.end(); it++){
		hAmino.push_back(it->first);
		//out << "amino acid:" << it->first << "; contain ptm: ";
		hNumPTM.push_back(it->second.size());
		for (int j = 0; j < it->second.size(); j++){
			//out << it->second[j] << " ";
			hPTM.push_back(it->second[j]);
		}
		//out << endl;
	}
	dPTM = hPTM;
	dNumPTM = hNumPTM;
	dAmino = hAmino;
	int aminoSize = hAmino.size();
	int n = x*y;
	int deltaBlocks = (n + 512 - 1) / 512;
	int numblocks = 32;
	int numthreads = n / numblocks + 1;

	cudaEvent_t start, stop;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	cudaEventRecord(start, NULL);
	//(int x, int y, vector<int>& diff, vector<int>& coords, vector<int>& numDeltaCoords, vector<int>& deltaCoordsx, vector<int>& deltaCoordsy, int tolerance, vector<int>& delta, int deltaSize)
	//CPU
	//cpuKernel2(x, y, diff, coords,numDeltaCoords,coordsx,coordsy,20000*tolerance,deltaInt,deltaSize);
	//GPU
	//t_KernelInitialDeltaDiag2 << <numblocks, numthreads >> >(x, y, thrust::raw_pointer_cast(&ddiff[0]), thrust::raw_pointer_cast(&dCoords[0]), thrust::raw_pointer_cast(&dNumDeltaCoords[0]), thrust::raw_pointer_cast(&dDeltaCoordsx[0]), thrust::raw_pointer_cast(&dDeltaCoordsy[0]), 20000 * tolerance, thrust::raw_pointer_cast(&ddeltaInt[0]), deltaSize);
	//KernelInitialDeltaDiagInt<<<numblocks,numthreads>>>(x, y, thrust::raw_pointer_cast(&dmassInt[0]),thrust::raw_pointer_cast(&dpeaksInt[0]),thrust::raw_pointer_cast(&dCoords[0]),thrust::raw_pointer_cast(&dNumDeltaCoords[0]),thrust::raw_pointer_cast(&dDeltaCoordsx[0]),thrust::raw_pointer_cast(&dDeltaCoordsy[0]),20000*tolerance,thrust::raw_pointer_cast(&ddeltaInt[0]),deltaSize);

	//KernelInitialDeltaDiag<<<numblocks,numthreads>>>(x, y, thrust::raw_pointer_cast(&dmass[0]),thrust::raw_pointer_cast(&dpeaks[0]),thrust::raw_pointer_cast(&dCoords[0]),thrust::raw_pointer_cast(&dNumDeltaCoords[0]),thrust::raw_pointer_cast(&dDeltaCoordsx[0]),thrust::raw_pointer_cast(&dDeltaCoordsy[0]),2*tolerance,thrust::raw_pointer_cast(&dDelta[0]),deltaSize);
	cudaEventRecord(stop, NULL);
	cudaEventSynchronize(start);
	cudaEventSynchronize(stop);

	float msecTotal = 0.0f;
	cudaEventElapsedTime(&msecTotal, start, stop);
	//std::cout << "diag_delta并行时间:" << msecTotal << "ms" << endl;

	hNumDeltaCoords = dNumDeltaCoords;
	hDeltaCoordsx = dDeltaCoordsx;
	hDeltaCoordsy = dDeltaCoordsy;
	*/

	//out << "--------------------------------------------" << endl;
	//out << "--------------------------------------------" << endl;
	/*
	int testNumdelat = 0;
	for (int i = 0; i < hNumDeltaCoords.size(); i++){
		if (hNumDeltaCoords[i] != 0){
			int tempi = i%x;
			int tempj = i / x;
			testNumdelat += hNumDeltaCoords[i];
			//out << "(" << cutMasses[tempi] << "," << cutMonoPeaks[tempj] << ") ";
			//out << "(" << tempi << "," << tempj << ") -> ";
			for (int k = 0; k < 5; k++){
				int tx = hDeltaCoordsx[5 * i + k];
				int ty = hDeltaCoordsy[5 * i + k];
				//out << "(" << cutMasses[tx] << "," << cutMonoPeaks[ty] << ") ";
				//out << "(" << tx << "," << ty << "), ";
			}
			//out << endl;
		}
	}
	*/
	
}