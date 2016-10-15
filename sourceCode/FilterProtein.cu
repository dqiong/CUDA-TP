#include "FilterProtein.cuh"

double findKey(double b, map<double, int>& m){
	double result = -1;
	map<double, int>::iterator it;
	for (it = m.begin(); it!= m.end(); it++){
		if (fabs(b - it->first) < 0.015){
			result = it->second;
			break;
		}
	}
	return result;

}
void calWeight(const vector<double>& peaks, const vector<double>& masses){
	clock_t start = clock();
	map<double, int> crossLine;
	for (int i = 0; i < masses.size(); i++){
		for (int j = 0; j < peaks.size(); j++){
			double b = peaks[j] - masses[i];
			double d = findKey(b, crossLine);
			if (d==-1){
				crossLine.insert(pair<double, int>(b, 1));
			}
			else{
				crossLine[d]++;
			}
		}
	}
	clock_t end = clock();
	cout << "一次计算分数时间(自己的版本):" << end - start << "ms" << endl;
}

void approxConvolution(const vector<double>& a, const vector<double>& b, double e){
	clock_t start = clock();
	int n = a.size();
	int m = b.size();
	double ma = a[n - 1];
	double mb = b[m - 1];
	int max = (ma + mb) / e;
	int score = 0;
	vector<int> c(max+2, 0);
	for (int i = 0; i < n; i++){
		for (int j = 0; j < m; j++){
			int d = (b[j] - a[i]+ma) / e;
			c[d] = c[d] + 1;
		}
	}
	for (int i = 0; i < max; i++){
		int k = c[i] + c[i + 1];
		if (k>score)
			score = k;
	}
	clock_t end = clock();
	//cout << "一次计算分数时间(approx版本1):" << end - start << "ms" << endl;
	//cout << "score:" << score << endl;
}
int approxConvolution2(const vector<double>& masses, const vector<double>& peaks, double e){
	//clock_t start = clock();
	int n = masses.size();
	int m = peaks.size();
	double mb = peaks[m - 1];
	double ma = masses[n - 1];
	int score = 0;
	int max = (ma + mb) / e;
	//vector<int> c(max+2, 0);
	
	int *c = new int[max+2];
	memset(c,0,(max+2)*sizeof(int));

	for (int i = 0; i < n; i++){
		for (int j = 0; j < m; j++){
			int d = (peaks[j] - masses[i]+ma) / e;
			//cout << "d::" << d <<" m-j:"<<m-j<< endl;
			c[d] = c[d] + 1;
			if (c[d]>score){
				score = c[d];
			}
			d++;
			c[d] = c[d] + 1;
			if (c[d]>score){
				score = c[d];
			}
		}
	}
	free(c);
	return score;
	//clock_t end = clock();
	//std::cout << "一次计算分数时间(approx版本2):" << end - start << "ms" << endl;
	//cout << "p值:" << dmax << endl;
	//std::cout << "method2 score:" << score << endl;

}
vector<int> restrictedConvolution(vector<double>& a, vector<double>& b, double p, double delta, double e){
	int n = a.size();
	int m = b.size();
	double ma = a[n - 1];
	double mb = b[m - 1];
	int max = (ma + mb) / e;
	vector<int> c(max+2, 0);
	int i = 0, j = 0;
	while (i < n&&j < m){
		int d = b[j] - a[i];
		if (d < p){
			j++;
		}
		else{
			i++;
			if (d < p + delta)
				c[(d + ma) / e]++;
		}
	}
	return c;
}
void exactConvolution(const vector<double>& a, const vector<double>& b, double e){
	clock_t c_start = clock();
	int n = a.size();
	int m = b.size();
	vector<double> c(m*n, 0);
	for (int i = 0; i < n; i++){         //calculate the convolution array at resolution ee>e
		for (int j = 0; j < m; j++){
			c.push_back(b[j] - a[i]);
		}
	}
	sort(c.begin(), c.end());
	//cout << c[c.size()-1] << endl;
	int score = 0;
	int k = 0;
	for (int i = 0; i < c.size(); i++){
		k = 0;
		double p = c[i];
		if (p != 0){
			for (int j = i + 1; j < c.size() && (j - i) < 50; j++){
				if (c[j] < (p + e)){
					k++;
					//cout << c[i + j] << " " << c[i] + e << endl;
				}

			}
			if (k>score){
				score = k;
				//cout << p << "  " << score << endl;
			}
		}
		
			
	}
	clock_t c_end = clock();
	//cout << "一次计算分数时间(exact:):" << c_end - c_start << "ms" << endl;
	//cout << "exact score:" << score << endl;
}
int msFilter(const vector<double>& a, const vector<double>& b, double e){
	//first convolution stage
	double amino_acid_min = 75.06;
	//clock_t c_start = clock();
	int t = 5;
	int n = a.size();
	int m = b.size();
	double ma = a[n - 1];   // maximum protein  weight
	double mb = b[m - 1];   //maximum spectrum
	double ee = 1;
	int maxee = (ma + mb) ;
	vector<double> begin;
	vector<double> end;
	vector<double> diff;
	//vector<int> index;
	vector<int> tempd;
	//vector<int> cc(maxee+2, 0);
	int *cc = new int[maxee+2];
	memset(cc,0,(maxee+2)*sizeof(int));
	for (int i = 0; i < n; i++){         //calculate the convolution array at resolution ee>e
		double protmass = a[i] + mb;
		for (int j = 0; j < m; j++){
			int d = protmass-b[j];
			cc[d]++;
			if (cc[d]>t){
				//t = cc[d];
				if (find(tempd.begin(), tempd.end(), d) == tempd.end()){
					//index.push_back(j*n + i);
					tempd.push_back(d);
					diff.push_back( a[i]-b[j]);
				}		
			}
		}
	}
	//cout << "1Da num:" << index.size() << endl;

	//stage2
	//clock_t stage2s = clock();
	sort(diff.begin(), diff.end());
	
	for (int i = 0; i < diff.size(); i++){
		if (i+1<diff.size()&&diff[i + 1] - diff[i]+2*ee < amino_acid_min){
			begin.push_back(diff[i] - ee);
			end.push_back(diff[i + 1] + ee);
			i++;
		}
		else{
			begin.push_back(diff[i] - ee);
			end.push_back(diff[i] + ee);
		}
	}
	
	t = 0;
	//int max = (ma + mb) / e;
	//vector<int> c(max+2, 0);
	//int *c = new int[max+2];
	//memset(c,0,(max+2)*sizeof(int));
	//clock_t stage2e = clock();
	for (int k = 0; k < begin.size(); k++){

		int i = 0, j = 0;
		double first = begin[k];
		double last = end[k];
		int limit = ((last - first) / e) + 2;
		int *c = new int[limit];
		memset(c, 0, (limit)*sizeof(int));
		/*
		if (k>0 && diff[index[k]] - diff[index[k - 1]] < ee){
		p = diff[index[k - 1]];
		delta = diff[index[k]] - diff[index[k - 1]];
		}*/
		//cout << diff[index[k]] << " ";

		while (i < n&&j < m){
			double d = a[i] - b[j];
			if (d >= last){
				j++;
			}
			else{
				i++;

				if (d>first){
					//  here important
					int indexd = (d - first) / e;
					c[indexd] ++;
					/*
					if (c[indexd]>t){
					t = c[indexd];
					}
					indexd++;
					c[indexd] ++;
					if (c[indexd] > t){
					t = c[indexd];
					}
					*/
				}
			}
		}
		int sum = c[limit - 1]+c[limit-2];
		int pos = limit - 3;
		while (pos >= 0){
			sum += c[pos];
			sum -= c[pos + 2];
			pos--;
			if (sum > t)
				t = sum + 1;
		}
		free(c);
	}
	free(cc);
	return t;
	//return t;
	//clock_t c_end = clock();
	//cout << "stage2:" << stage2e - stage2s << endl;
	//cout << "一次计算分数时间(msFilter:):" << c_end - c_start << "ms" << endl;
	//cout << "msFilter score:" << t << endl;
	//cout << "score:" << score << endl;

}
void msFilterOneSpectra(const vector<vector<double>>& a, const vector<double>& b, double e){
	//first convolution stage
	double amino_acid_min = 75.06;
	//clock_t c_start = clock();
	for (int idx_mass = 0; idx_mass < a.size(); idx_mass++){
		int t = 5;
		int n = a[idx_mass].size();
		int m = b.size();
		double ma = a[idx_mass][n - 1];
		double mb = b[m - 1];
		double ee = 1;
		int maxee = (ma + mb) / ee;
		vector<double> begin;
		vector<double> end;
		vector<double> diff;
		//vector<int> index;
		vector<int> tempd;
		//vector<int> cc(maxee+2, 0);
		int *cc = new int[maxee + 2];
		memset(cc, 0, (maxee + 2)*sizeof(int));
		for (int i = 0; i < n; i++){         //calculate the convolution array at resolution ee>e
			for (int j = 0; j < m; j++){
				int d = (b[j] - a[idx_mass][i] + ma) / ee;
				cc[d] = cc[d] + 1;
				if (cc[d]>t){
					t = cc[d];
					if (find(tempd.begin(), tempd.end(), d) == tempd.end()){
						//index.push_back(j*n + i);
						tempd.push_back(d);
						diff.push_back(b[j] - a[idx_mass][i]);
					}
				}
			}
		}
		//cout << "1Da num:" << index.size() << endl;

		//stage2
		clock_t stage2s = clock();
		sort(diff.begin(), diff.end());

		for (int i = 0; i < diff.size(); i++){
			if (i + 1 < diff.size() && diff[i + 1] - diff[i] + 2 * ee < amino_acid_min){
				begin.push_back(diff[i] - ee);
				end.push_back(diff[i + 1] + ee);
				i++;
			}
			else{
				begin.push_back(diff[i] - ee);
				end.push_back(diff[i] + ee);
			}
		}

		t = 0;
		int max = (ma + mb) / e;
		//vector<int> c(max+2, 0);
		int *c = new int[max + 2];
		memset(c, 0, (max + 2)*sizeof(int));
		//clock_t stage2e = clock();
		for (int k = 0; k < begin.size(); k++){

			int i = 0, j = 0;
			double first = begin[k];
			double last = end[k];
			/*
			if (k>0 && diff[index[k]] - diff[index[k - 1]] < ee){
			p = diff[index[k - 1]];
			delta = diff[index[k]] - diff[index[k - 1]];
			}*/
			//cout << diff[index[k]] << " ";

			while (i < n&&j < m){
				double d = b[j] - a[idx_mass][i];
				if (d < first){
					j++;
				}
				else{
					i++;
					if (d < last){
						//  cout << d << endl;
						int indexd = (d + ma) / e;
						c[indexd] = c[indexd] + 1;;
						if (c[indexd]>t){
							t = c[indexd];
						}

						indexd++;
						c[indexd] = c[indexd] + 1;
						if (c[indexd] > t){
							t = c[indexd];
						}

					}

				}
			}
		}
		free(cc);
		free(c);
		//return t;
		//clock_t c_end = clock();
		//cout << "stage2:" << stage2e - stage2s << endl;
		//cout << "一次计算分数时间(msFilter stage1:):" << c_end - c_start << "ms" << endl;
		//cout << "cpu msFilter:" << t << endl;
		//cout << "score:" << score << endl;
	}
	

}
/*
void filterProtein(vector<vector<mass_t>>& cutMasses,vector<MonoSpectrum>& spectrum,vector<vector<mass_t>>& cutPeaksMonoMass){
	clock_t msfilter_start = clock();
	char *path = "output\\time_result.txt";
	ofstream out(path);
	int limit = 800;
	//double parentMass2 = spectrum[2].getParentMonoMass();
	int spectraSize = min(spectrum.size(), cutPeaksMonoMass.size());
	vector<double> parentMasses;
	for (int i = 0; i < spectraSize ; i++){
		parentMasses.push_back(spectrum[i].getParentMonoMass());
	}
	vector<vector<int>> index(spectraSize );
	int filterAllSize = 0;
	for (int i = 0; i < spectraSize ; i++){
		for (int j = 0; j < cutMasses.size(); j++){
			if (fabs(parentMasses[i] - cutMasses[j][cutMasses[j].size() - 1]) < limit)
				index[i].push_back(j);
		}
		filterAllSize += index[i].size();
	}
	out << "每个谱的平均蛋白质个数:" << filterAllSize / spectraSize << endl;
	
	for (int i = 0; i<spectraSize ; i++){
		for (int j = 0; j < index[i].size(); j++){
			approxConvolution2(cutMasses[index[i][j]], cutPeaksMonoMass[i], 0.015);
		}
	}

	Config			config;
	config.init_with_defaults();
	vector<PTM> list = config.getAllPTMs();
	vector<mass_t> ptm;
	for (int j = 0; j < list.size(); j++){
		ptm.push_back(list[j].delta);
	}
	thrust::device_vector<double> ddelta = ptm;


	clock_t msfilter_end = clock();
	out << "ms-filter time:" << msfilter_end - msfilter_start << "ms" << endl;
	out<<"ms-filter time:" << (msfilter_end - msfilter_start)/60000 << "min" << endl;
	//cout << "过滤出:"<<index.size() << endl;

}
*/
void testTime(vector<vector<mass_t>>& cutMasses, vector<vector<mass_t>>& cutPeaksMonoMass){
	int proteinSize = cutMasses.size();
	double e = 0.2;
	clock_t begin_exact = clock();
	for (int i = 0; i < proteinSize; i++){
		exactConvolution(cutMasses[i], cutPeaksMonoMass[5],e);
	}
	clock_t end_exact = clock();
	cout << "exact time:" << end_exact - begin_exact << "ms" << endl;

	clock_t begin_app = clock();
	for (int i = 0; i < proteinSize; i++){
		approxConvolution(cutMasses[i], cutPeaksMonoMass[5],e);
	}
	clock_t end_app = clock();
	cout << "approx time:" << end_app - begin_app<< "ms" << endl;

	clock_t begin_app2 = clock();
	for (int i = 0; i < proteinSize; i++){
		approxConvolution2(cutMasses[i], cutPeaksMonoMass[5],e);
	}
	clock_t end_app2 = clock();
	cout << "approx2 time:" << end_app2 - begin_app2<< "ms" << endl;

	clock_t begin_ms = clock();
	for (int i = 0; i < proteinSize; i++){
		msFilter(cutMasses[i], cutPeaksMonoMass[5],e);
	}
	clock_t end_ms = clock();
	cout << "ms-filter time:" << end_ms - begin_ms<< "ms" << endl;

}
/*
void gpuFilter(vector<vector<mass_t>>& cutMasses,vector<vector<mass_t>>& cutPeaksMonoMass){
	thrust::host_vector<double> host_cutPeaks0;
	for (int i = 0; i < cutPeaksMonoMass[0].size(); i++){
		host_cutPeaks0.push_back(cutPeaksMonoMass[0][i]);
	}
	thrust::host_vector<int> cutMass_begin;
	thrust::host_vector<int> cutMass_end;
	thrust::host_vector<double> cutMassesTemp;
	thrust::device_vector<double> device_cutPeaks0 = host_cutPeaks0;
	int index = 0;
	for (int i = 0; i < cutMasses.size(); i++){
		cutMass_begin.push_back(index);
		for (int j = 0; j < cutMasses[i].size(); j++){
			cutMassesTemp.push_back(cutMasses[i][j]);
		}
		index += cutMasses[i].size() - 1;
		cutMass_end.push_back(index);
		index++;
	}
	
	thrust::device_vector<int> device_cutMass_begin(cutMass_begin.size());
	thrust::device_vector<int> device_cutMass_end(cutMass_end.size());
	thrust::device_vector<double> device_cutMasses(cutMassesTemp.size());
	device_cutMass_begin = cutMass_begin;
	device_cutMass_end = cutMass_end;
	device_cutMasses = cutMassesTemp;

	char *path = "output\\time_result.txt";
	ofstream out(path);

	//cout << cutMassesTemp[cutMass_end[6]] << "  " << host_cutPeaks0[host_cutPeaks0.size() - 1] << endl;
	//cout << "5 test:"<<(cutMassesTemp[cutMass_end[6]]+ host_cutPeaks0[host_cutPeaks0.size()-1])/0.015<< endl;
	out << "device_cutMass_begin.size() " << device_cutMass_begin.size() << endl;
	out << "device_cutMass_end.size() " << device_cutMass_end.size() << endl;
	out << "device_cutMasses.size() " << device_cutMasses.size() << endl;
	out << "device_cutPeaks0.size() " << device_cutPeaks0.size() << endl;

	clock_t dev_start = clock();
	//cudaEvent_t device_start0, device_stop0;
	//cudaEventCreate(&device_start0);
	//cudaEventCreate(&device_stop0);
	//cudaEventRecord(device_start0, NULL);
	//cout << "here" << endl;
	//testt << <1, 1 >> >();
	//test protein filter GPU
	//thrust::device_vector<double> device_cutMass2 = cutMasses[2];

	cudaDeviceSetLimit(cudaLimitMallocHeapSize, 100 * 1024*1024);
	//gpuFilter << <1, 3 >> >();
	//gpuFilter << <1,  1 >> >(thrust::raw_pointer_cast(&device_cutMass_begin[0]), thrust::raw_pointer_cast(&device_cutMass_end[0]), thrust::raw_pointer_cast(&device_cutMasses[0]), thrust::raw_pointer_cast(&device_cutPeaks0[0]), device_cutPeaks0.size(),proteins.size());
	//gpuApproxConvolution<<<1,1>>>(thrust::raw_pointer_cast(&device_cutMass_begin[0]), thrust::raw_pointer_cast(&device_cutMass_end[0]), thrust::raw_pointer_cast(&device_cutMasses[0]), thrust::raw_pointer_cast(&device_cutPeaks0[0]), device_cutPeaks0.size(),proteins.size());
	//cudaSetDeviceFlags(cudaDeviceMapHost);
	//const int fi = 1;
	//cudaStream_t stream[fi];
	//for (int i = 0; i < fi; i++){
	//	cudaStreamCreate(&stream[i]);
		for (int j = 0; j < 10; j++){
			int testIndex = j;
			int ma_size = cutMasses[testIndex].size();
			int mb_size = cutPeaksMonoMass[0].size();
			thrust::device_vector<double> cutMassTest(ma_size);
			cutMassTest = cutMasses[testIndex];
			thrust::device_vector<double> cutPeakTest(mb_size);
			cutPeakTest = cutPeaksMonoMass[0];
			double ma = cutMasses[testIndex][ma_size - 1];
			double mb = cutPeaksMonoMass[0][mb_size - 1];
			int num = ma_size*mb_size;
			int max = ((ma + mb) / 0.015) + 2;
			int *c;
			//CUDA_CALL(cudaHostAlloc((void**)&c, max*sizeof(int), cudaHostAllocWriteCombined | cudaHostAllocMapped));
			cudaMalloc((void**)&c, sizeof(int)* max);
			int *dev_score;
			cudaMalloc((void**)&dev_score, sizeof(int)* 1);
			//int *dev_c;
			int score[1] = { 0 };
		//	int *hc=new int[max];
			cudaMemcpy(dev_score, score, sizeof(int), cudaMemcpyHostToDevice);
			//CUDA_CALL(cudaHostGetDevicePointer(&dev_c, c, 0));
			//cudaMemcpy(gpudata, data, sizeof(int) * DATA_SIZE,cudaMemcpyHostToDevice);
			calScore << <num / 256+1, 256 >> >(thrust::raw_pointer_cast(&cutMassTest[0]), thrust::raw_pointer_cast(&cutPeakTest[0]), ma_size, mb_size, ma, mb, 0.015, c, num, dev_score, max);
			//cudaThreadSynchronize();
			cudaMemcpy(score, dev_score, sizeof(int), cudaMemcpyDeviceToHost);
			//cudaMemcpy(hc, c, sizeof(int)*max, cudaMemcpyDeviceToHost);
			//int *tempresult = max_element(hc, hc + max);
			cout << "并行:" << score[0] << endl;
			//cout << "标记:" << testIndex << endl;
			cudaFreeHost(c);
			cudaFreeHost(dev_score);
		}
	//}
	//for (int i = 0; i < fi;i++)
	//	cudaStreamDestroy(stream[i]);
	clock_t dev_end = clock();
	//cudaEventRecord(device_stop0, NULL);
	//cudaEventSynchronize(device_start0);
	//cudaEventSynchronize(device_stop0);
	//float msecTotal0 = 0.0f;
	//cudaEventElapsedTime(&msecTotal0, device_start0, device_stop0);
	out << " Time GPU protein filter:" << dev_end-dev_start << "ms"<<endl;
}
__global__ void msFiltergpu(int* cutMass_begin, int* cutMass_end, double* a,double* b,int m){
	int idx = blockDim.x*blockIdx.x + threadIdx.x;
	//first convolution stage
	double amino_acid_min = 75.06;
	//clock_t c_start = clock();
	int t = 8;
	int index_start = cutMass_begin[idx];
	int index_end = cutMass_end[idx];
	int n = index_end-index_start+1;
	double ma = a[index_end];
	double mb = b[m - 1];
	double ee = 1;
	int maxee = (ma + mb) / ee;
	thrust::device_vector<double> begin;
	thrust::device_vector<double> end;
	thrust::device_vector<double> diff;
	//vector<int> index;
	thrust::device_vector<int> tempd;
	thrust::device_vector<int> cc(maxee+2, 0);

	for (int i = index_start; i < index_end+1; i++){         //calculate the convolution array at resolution ee>e
		for (int j = 0; j < m; j++){
			int d = (b[j] - a[i]+ma) / ee;
			cc[d] = cc[d] + 1;
			if (cc[d]>t){
				//find function
				if (!findVal(thrust::raw_pointer_cast(&tempd[0]),d,tempd.size())){
					//index.push_back(j*n + i);
					tempd.push_back(d);
					diff.push_back( b[j] - a[i]);
				}		
			}
		}
	}
	//cout << "1Da num:" << index.size() << endl;

	//stage2
	thrust::sort(diff.begin(), diff.end());
	double e = 0.015;
	for (int i = 0; i < diff.size(); i++){
		if (i+1<diff.size()&&diff[i + 1] - diff[i]+2*ee < amino_acid_min){
			begin.push_back(diff[i] - ee);
			end.push_back(diff[i + 1] + ee);
			i++;
		}
		else{
			begin.push_back(diff[i] - ee);
			end.push_back(diff[i] + ee);
		}
	}
	
	t = 0;
	int max = (ma + mb) / e;
	thrust::device_vector<int> c(max+2, 0);
	for (int k= 0; k < begin.size(); k++){
		
		int i = 0, j = 0;
		double first = begin[k];
		double last = end[k];
		//cout << diff[index[k]] << " ";
		
		while (i < n&&j < m){
		   double d = b[j] - a[index_start+i];
		   if (d < first){
			  j++;
		   }
		   else{
			   i++;
			   if (d < last){
				 //  cout << d << endl;
				   int indexd = (d + ma) / e;
				   c[indexd] = c[indexd] + 1;;
				   if (c[indexd]>t){
					   t = c[indexd];
				   }
				   
				  indexd++;
				   c[indexd]=c[indexd]+1;
				   if (c[indexd] > t){
					   t = c[indexd];
				   }
				   
			   }

		   }
		}
	}
	printf("gpu msFilter: %d", t);
	
}

__global__ void testt(){
	int idx = blockDim.x*blockIdx.x + threadIdx.x;
	printf("idx: %d", idx);
}
*/