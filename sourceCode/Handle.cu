#include "Handle.cuh"


//when there is no mass shift,the proceeds of array D and M
void InitialArray(thrust::host_vector<int>& d, thrust::host_vector<int>& m, int x, int y, thrust::host_vector<CutPeak> peaks, thrust::host_vector<mass_t> masses, double tolerance){

	for (int i = 0; i < x; i++){
		m[i] = 0;
		/*
		if (fabs(peaks[i].monoMass - masses[i]) < tolerance){

		}
		*/
		d[i] = 0;
	}
	for (int j = 0; j < y; j++){
		m[j*x] = 0;
		d[j*x] = 0;
	}
	int count = 0;
	for (int j = 1; j < y; j++)
	for (int i = 1; i < x; i++){
		int index = j*x + i;
		int upIndex = (j - 1)*x + i;
		int leftIndex = j*x + i - 1;
		if (fabs(peaks[j].monoMass - masses[i]) < tolerance){
			count++;
			d[index] = count;
		}
		else{
			d[index] = 0;
		}
		m[index] = max(d[index], max(m[upIndex], m[leftIndex]));
	}
}
void InitialArray(vector<int>& d, vector<int>& m, int x, int y, vector<CutPeak> peaks, vector<mass_t> masses, double tolerance){
	for (int i = 0; i < x; i++){
		m[i] = 0;
		/*
		if (fabs(peaks[i].monoMass - masses[i]) < tolerance){

		}
		*/
		d[i] = 0;
	}
	for (int j = 0; j < y; j++){
		m[j*x] = 0;
		d[j*x] = 0;
	}
	int count = 0;
	for (int j = 1; j < y; j++)
	for (int i = 1; i < x; i++){
		int index = j*x + i;
		int upIndex = (j - 1)*x + i;
		int leftIndex = j*x + i - 1;
		if (fabs(peaks[j].monoMass - masses[i]) < tolerance){
			count++;
			d[index] = count;
		}
		else{
			d[index] = 0;
		}
		m[index] = max(d[index], max(m[upIndex], m[leftIndex]));
	}
}
void InitialPreviousPoint(Node* &root, thrust::host_vector<mass_t>& A, thrust::host_vector<mass_t>& B, int x, int y, Point* match){

	for (int i = 0; i < x; i++)
	for (int j = 0; j < y; j++){
		double b = B[j] - A[i];
		Node* tmp = root->Find(root, b);
		if (tmp == NULL){
			root->Insert(root, b, i, j);
			match[j*x + i].x = 0;       // 0 denote there is no previous point
			match[j*x + i].y = 0;
		}
		else{
			match[j*x + i].x = tmp->xx;
			match[j*x + i].y = tmp->yy;
			tmp->xx = i, tmp->yy = j;
		}
		// cout << "there is ok:" <<"noMasses:"<<i<<"  noPeaks:"<<j<< endl;
	}
}
void InitialPreviousPoint(Node* &root, vector<mass_t>& A, vector<mass_t>& B, int x, int y, Point *match){
	for (int i = 0; i < x; i++)
	for (int j = 0; j < y; j++){
		double b = B[j] - A[i];
		Node* tmp = root->Find(root, b);
		if (tmp == NULL){
			root->Insert(root, b, i, j);
			match[j*x + i].x = 0;       // 0 denote there is no previous point
			match[j*x + i].y = 0;
		}
		else{
			match[j*x + i].x = tmp->xx;
			match[j*x + i].y = tmp->yy;
			tmp->xx = i, tmp->yy = j;
		}
		// cout << "there is ok:" <<"noMasses:"<<i<<"  noPeaks:"<<j<< endl;
	}
}
void SplitString(const string& s, vector<std::string>& v, const std::string& c)
{
	string::size_type pos1, pos2;
	pos2 = s.find(c);
	pos1 = 0;
	while (std::string::npos != pos2)
	{
		v.push_back(s.substr(pos1, pos2 - pos1));

		pos1 = pos2 + c.size();
		pos2 = s.find(c, pos1);
	}
	if (pos1 != s.length())
		v.push_back(s.substr(pos1));
}
void InitialMapAminoPtm(map<char, bool>& m){
	m.insert(make_pair('A', false));
	m.insert(make_pair('R', false));
	m.insert(make_pair('D', false));
	m.insert(make_pair('C', false));
	m.insert(make_pair('Q', false));
	m.insert(make_pair('E', false));
	m.insert(make_pair('H', false));
	m.insert(make_pair('I', false));
	m.insert(make_pair('G', false));
	m.insert(make_pair('N', false));
	m.insert(make_pair('L', false));
	m.insert(make_pair('K', false));
	m.insert(make_pair('M', false));
	m.insert(make_pair('F', false));
	m.insert(make_pair('P', false));
	m.insert(make_pair('S', false));
	m.insert(make_pair('T', false));
	m.insert(make_pair('W', false));
	m.insert(make_pair('Y', false));
	m.insert(make_pair('V', false));
}

void gpuFunction(vector<vector<mass_t>>& cutMasses,vector<CutPeak>& cutPeaks,vector<mass_t>& cutPeaksMonoMass,vector<mass_t> ptm, vector<int>& tscore){
	mass_t ppm = 15.0;
	int total = 10<cutMasses.size()?10:cutMasses.size();
	//const int total = temp_total;
	unsigned int numSpecific = 6;
	double tolerance = ppm / 1000;

	Config			config;
	
	/*
	char *path1 = "output\\result.txt";
	ofstream out(path1);
	out.setf(ios::showpoint);
	out.precision(6);
	out.setf(ios::fixed);
	*/

	clock_t start = clock();
	thrust::device_vector<double> ddelta = ptm;
	
	int numPeaks=cutPeaksMonoMass.size();
	int *numMasses=new int[total];
	int *n=new int[total];
	vector<vector<int>> v_D0(total);
	vector<vector<int>> v_M0(total);
	clock_t c1 = clock();
	for (int i = 0; i < total; i++){
		numMasses[i] = cutMasses[i].size();
		n[i] = numMasses[i]*numPeaks;
		v_D0[i].resize(n[i]);
		v_M0[i].resize(n[i]);
		InitialArray(v_D0[i], v_M0[i], numMasses[i], numPeaks, cutPeaks, cutMasses[i], tolerance);
	}

	vector<thrust::device_vector<int>> device_D0(total);
	vector<thrust::device_vector<int>> device_M0(total);
	vector<thrust::device_vector<int>> device_D1(total);
	vector<thrust::device_vector<int>> device_M1(total);
	//PreviousDiag
	vector<thrust::device_vector<int>> dpreviousMatchx(total);
	vector<thrust::device_vector<int>> dpreviousMatchy(total);
	//PreviousDeltaDiag
	vector<thrust::device_vector<int>> dNumDeltaCoords(total);
	vector<thrust::device_vector<int>> dDeltaCoordsx(total);
	vector<thrust::device_vector<int>> dDeltaCoordsy(total);
	//stream
	cudaStream_t stream[10];
	for (int i = 0; i < total; i++){

		cudaStreamCreate(&stream[i]);

		device_D0[i]=v_D0[i];
		device_M0[i]=v_M0[i];
	//	cudaMemcpyAsync(thrust::raw_pointer_cast(&device_D0[i][0]), &v_D0[i][0], n[i]*sizeof(int), cudaMemcpyHostToDevice, stream[i]);
		//cudaMemcpyAsync(thrust::raw_pointer_cast(&device_M0[i][0]), &v_M0[i][0], n[i]*sizeof(int), cudaMemcpyHostToDevice, stream[i]);

		Point* previousMatch = new Point[n[i]];
		for (int j = 0; j < numPeaks; j++)
		for (int k = 0; k < numMasses[i]; k++){
			previousMatch[j*numMasses[i] + k].x = 0;
			previousMatch[j*numMasses[i] + k].y = 0;
		}
		Node* root;
		root = NULL;
		InitialPreviousPoint(root, cutMasses[i], cutPeaksMonoMass, numMasses[i], numPeaks, previousMatch);

		//PreviousDiag
		dpreviousMatchx[i].resize(n[i]);
		dpreviousMatchy[i].resize(n[i]);
		//PreviousDeltaDiag
		dNumDeltaCoords[i].resize(n[i]);
		dDeltaCoordsx[i].resize(5 * n[i], 0);
		dDeltaCoordsy[i].resize(5 * n[i], 0);


		for (int j = 0; j < numPeaks; j++){
			for (int k = 0; k < numMasses[i]; k++){
				int previous_x = previousMatch[j*numMasses[i] + k].x;
				int previous_y = previousMatch[j*numMasses[i] + k].y;
				dpreviousMatchx[i][j*numMasses[i] + k] = previous_x;
				dpreviousMatchy[i][j*numMasses[i] + k] = previous_y;


			}
		}

		device_D1[i].resize(n[i], 0);
		device_M1[i].resize(n[i], 0);
		root->Delete(root);
	}

	//cout << "All right herhe" << endl;

	vector<vector<int>> diff(total);
	vector<vector<int>> massInt(total);
	vector<vector<int>> peakInt(total);
	vector<vector<int>> deltaInt(total);
	vector<thrust::device_vector<int>> d_diff(total);
	vector < thrust::device_vector<int>> ddeltaInt(total);

	for (int i = 0; i < total; i++){

		int numblocks = 32;
	    int numdeltathreads = n[i] / numblocks + 1;

		diff[i].resize(n[i]);
		convert(cutMasses[i], massInt[i], cutPeaksMonoMass, peakInt[i], ptm, deltaInt[i]);
		calcDiff(massInt[i], peakInt[i], diff[i]);

		d_diff[i]=diff[i];
		ddeltaInt[i]=deltaInt[i];
		//cudaMemcpyAsync(thrust::raw_pointer_cast(&d_diff[i][0]), &diff[i][0], diff[i].size()*sizeof(int), cudaMemcpyHostToDevice, stream[i]);
		//cudaMemcpyAsync(thrust::raw_pointer_cast(&ddeltaInt[i][0]), &deltaInt[i][0], deltaInt[i].size()*sizeof(int), cudaMemcpyHostToDevice, stream[i]);


	/*	cudaEvent_t device_start0, device_stop0;
		cudaEventCreate(&device_start0);
		cudaEventCreate(&device_stop0);
		cudaEventRecord(device_start0, NULL);*/
		KernelInitialDeltaDiag2 << <numblocks, numdeltathreads ,0,stream[i]>> >(numMasses[i], numPeaks, thrust::raw_pointer_cast(&d_diff[i][0]), thrust::raw_pointer_cast(&dNumDeltaCoords[i][0]), thrust::raw_pointer_cast(&dDeltaCoordsx[i][0]), thrust::raw_pointer_cast(&dDeltaCoordsy[i][0]), 20000 * tolerance, thrust::raw_pointer_cast(&ddeltaInt[i][0]), deltaInt[i].size());
		/*cudaEventRecord(device_stop0, NULL);
		cudaEventSynchronize(device_start0);
		cudaEventSynchronize(device_stop0);
		float msecTotal0 = 0.0f;
		cudaEventElapsedTime(&msecTotal0, device_start0, device_stop0);
		out << " GPU计算delta_diag:" << msecTotal0 << "ms";*/
		//KernelInitialDeltaDiag << <numblocks, numdeltathreads >> >(numMasses, numPeaks, thrust::raw_pointer_cast(&d_diff[0]), thrust::raw_pointer_cast(&dNumDeltaCoords[0]), thrust::raw_pointer_cast(&dDeltaCoordsx[0]), thrust::raw_pointer_cast(&dDeltaCoordsy[0]), 2 * tolerance, thrust::raw_pointer_cast(&ddelta[0]), ddelta.size());

		int internalCellRow = 4;
		int internelCellColumn = 4;
		int blocks = 1, threads = numPeaks;
		/*cudaEvent_t device_start, device_stop;
		cudaEventCreate(&device_start);
		cudaEventCreate(&device_stop);
		cudaEventRecord(device_start, NULL);*/
		for (int k = 0; k < numSpecific / 2; k++){

			gpu_dp_kernel2 << <blocks, threads,0,stream[i] >> >(numMasses[i], numPeaks, thrust::raw_pointer_cast(&device_M0[i][0]), thrust::raw_pointer_cast(&device_D0[i][0]), thrust::raw_pointer_cast(&device_M1[i][0]), thrust::raw_pointer_cast(&device_D1[i][0]), thrust::raw_pointer_cast(&dpreviousMatchx[i][0]), thrust::raw_pointer_cast(&dpreviousMatchy[i][0]), thrust::raw_pointer_cast(&dNumDeltaCoords[i][0]), thrust::raw_pointer_cast(&dDeltaCoordsx[i][0]), thrust::raw_pointer_cast(&dDeltaCoordsy[i][0]), internalCellRow, internelCellColumn, numMasses[i] / internalCellRow, numMasses[i]%internalCellRow, numPeaks / internelCellColumn);

			gpu_dp_kernel2 << <blocks, threads,0,stream[i] >> >(numMasses[i], numPeaks, thrust::raw_pointer_cast(&device_M1[i][0]), thrust::raw_pointer_cast(&device_D1[i][0]), thrust::raw_pointer_cast(&device_M0[i][0]), thrust::raw_pointer_cast(&device_D0[i][0]), thrust::raw_pointer_cast(&dpreviousMatchx[i][0]), thrust::raw_pointer_cast(&dpreviousMatchy[i][0]), thrust::raw_pointer_cast(&dNumDeltaCoords[i][0]), thrust::raw_pointer_cast(&dDeltaCoordsx[i][0]), thrust::raw_pointer_cast(&dDeltaCoordsy[i][0]), internalCellRow, internelCellColumn, numMasses[i] / internalCellRow, numMasses[i]%internalCellRow, numPeaks/ internelCellColumn);

		}
		/*cudaThreadSynchronize();
		cudaEventRecord(device_stop, NULL);
		cudaEventSynchronize(device_start);
		cudaEventSynchronize(device_stop);
		float msecTotal = 0.0f;
		cudaEventElapsedTime(&msecTotal, device_start, device_stop);
		out << " GPU计算6次:" << msecTotal << "ms";
		out << "   no:" << i << ";GPU specific PTM 5: " << device_M1[n - 1] << endl;*/

	}
	cudaThreadSynchronize();
	clock_t end = clock();
	for (int i = 0; i < total; i++){
		cudaStreamDestroy(stream[i]);
		vector<int>(diff[i]).swap(diff[i]);
		vector<int>(massInt[i]).swap(massInt[i]);
		vector<int>(peakInt[i]).swap(peakInt[i]);
		vector<int>(deltaInt[i]).swap(deltaInt[i]);
		vector<int>(v_D0[i]).swap(v_D0[i]);
		vector<int>(v_M0[i]).swap(v_M0[i]);

		thrust::device_vector<int >(d_diff[i]).swap(d_diff[i]);
		thrust::device_vector<int >(ddeltaInt[i]).swap(ddeltaInt[i]);
		thrust::device_vector<int >(device_D0[i]).swap(device_D0[i]);
		thrust::device_vector<int >(device_M0[i]).swap(device_M0[i]);
		thrust::device_vector<int >(device_D1[i]).swap(device_D1[i]);
		thrust::device_vector<int >(device_M0[i]).swap(device_M0[i]);
		thrust::device_vector<int >(dpreviousMatchx[i]).swap(dpreviousMatchx[i]);
		thrust::device_vector<int >(dpreviousMatchy[i]).swap(dpreviousMatchy[i]);
		thrust::device_vector<int >(dNumDeltaCoords[i]).swap(dNumDeltaCoords[i]);
		thrust::device_vector<int >(dDeltaCoordsx[i]).swap(dDeltaCoordsx[i]);
		thrust::device_vector<int >(dDeltaCoordsy[i]).swap(dDeltaCoordsy[i]);
		//out << "   no:" << i << ";GPU specific PTM 6: " << device_M1[i][n[i] - 1] << endl;
		tscore.push_back(device_M1[i][n[i] - 1]);
	}

	//out << "--------------------------------------------" << endl;
	//out << "--------------------------------------------" << endl;
	//out << "Total time:" << end - start << "ms" << endl;
}

void gpuFunction2(){
	char peak_list_path[10][512];
	char sequence_path[512];
	unsigned int numSpecific = 6;
	unsigned int numGeneral = 0;
	unsigned int numForms = 1;
	mass_t ppm = 15.0;
	const int total = 10;

	double tolerance = ppm / 1000;

	Config			config;
	MonoSpectrum	spectrum[total];
	FragmentSet		fragments[total];
	ProteinSequence sequence;

	char *path1 = "output\\result.txt";
	ofstream out(path1);
	out.setf(ios::showpoint);
	out.precision(6);
	out.setf(ios::fixed);

	strcpy(sequence_path, "data\\Hela_seq.txt");
	config.init_with_defaults();
	vector<PTM> list = config.getAllPTMs();
	vector<mass_t> ptm;
	for (int j = 0; j < list.size(); j++){
		ptm.push_back(list[j].delta);
	}
	thrust::device_vector<double> ddelta = ptm;
	sequence.read(sequence_path, &config);
	vector<mass_t> cutMasses;
	sequence.calcNonModifiedCutMasses(cutMasses);

	string str = "";
	string str0 = "data\\HELA_spec";
	string str1 = "";
	string str2 = "_mono.txt";
	stringstream ss;

	clock_t start = clock();
	vector<vector<CutPeak>> cutPeaks(total);
	for (int i = 0; i < total; i++){
		ss << i + 1;
		str1 = ss.str();
		//ss.str("");
		str.append(str0).append(str1).append(str2);
		ss.str("");
		strcpy(peak_list_path[i], str.c_str());
		str = "";

		fragments[i].initECD();
		spectrum[i].readAllPeakList(peak_list_path[i], ppm, ppm);
		//从文件中读出峰值并转化为cutpeak
		spectrum[i].convertToCutPeaks(fragments[i], cutPeaks[i]);
	}

	vector<vector<mass_t>> cutPeaksMonoMass(total);
	for (int i = 0; i < cutPeaks.size(); i++){
		for (int j = 0; j < cutPeaks[i].size(); j++){
			cutPeaksMonoMass[i].push_back(cutPeaks[i][j].monoMass);
		}
	}
	//蛋白质过滤测试
	//calWeight(cutPeaksMonoMass[0], cutMasses);
	for (int i = 0; i < 10; i++){
		cout << "number." << i << " spectrum" << endl;
		exactConvolution(cutMasses, cutPeaksMonoMass[i], 0.015);
		approxConvolution(cutMasses, cutPeaksMonoMass[i], 0.015);
		approxConvolution2(cutMasses, cutPeaksMonoMass[i], 0.015);
		msFilter(cutMasses, cutPeaksMonoMass[i], 0.015);
		cout << "--------------------------------------" << endl;
	}
	
	//cout << "All right!!" << endl;

	int numPeaks[total];
	int numMasses = cutMasses.size();
	int n[total];
	vector<vector<int>> v_D0(total);
	vector<vector<int>> v_M0(total);
	clock_t c1 = clock();
	for (int i = 0; i < total; i++){
		numPeaks[i] = cutPeaks[i].size();
		n[i] = numMasses*numPeaks[i];
		v_D0[i].resize(n[i]);
		v_M0[i].resize(n[i]);
		InitialArray(v_D0[i], v_M0[i], numMasses, numPeaks[i], cutPeaks[i], cutMasses, tolerance);
	}



	vector<thrust::device_vector<int>> device_D0(total);
	vector<thrust::device_vector<int>> device_M0(total);
	vector<thrust::device_vector<int>> device_D1(total);
	vector<thrust::device_vector<int>> device_M1(total);
	//PreviousDiag
	vector<thrust::device_vector<int>> dpreviousMatchx(total);
	vector<thrust::device_vector<int>> dpreviousMatchy(total);
	//PreviousDeltaDiag
	vector<thrust::device_vector<int>> dNumDeltaCoords(total);
	vector<thrust::device_vector<int>> dDeltaCoordsx(total);
	vector<thrust::device_vector<int>> dDeltaCoordsy(total);
	//stream
	cudaStream_t stream[total];
	for (int i = 0; i < total; i++){

		cudaStreamCreate(&stream[i]);

		device_D0[i]=v_D0[i];
		device_M0[i]=v_M0[i];
	//	cudaMemcpyAsync(thrust::raw_pointer_cast(&device_D0[i][0]), &v_D0[i][0], n[i]*sizeof(int), cudaMemcpyHostToDevice, stream[i]);
		//cudaMemcpyAsync(thrust::raw_pointer_cast(&device_M0[i][0]), &v_M0[i][0], n[i]*sizeof(int), cudaMemcpyHostToDevice, stream[i]);

		Point* previousMatch = new Point[n[i]];
		for (int j = 0; j < numPeaks[i]; j++)
		for (int k = 0; k < numMasses; k++){
			previousMatch[j*numMasses + k].x = 0;
			previousMatch[j*numMasses + k].y = 0;
		}
		Node* root;
		root = NULL;
		InitialPreviousPoint(root, cutMasses, cutPeaksMonoMass[i], numMasses, numPeaks[i], previousMatch);
		//PreviousDiag
		dpreviousMatchx[i].resize(n[i]);
		dpreviousMatchy[i].resize(n[i]);
		//PreviousDeltaDiag
		dNumDeltaCoords[i].resize(n[i]);
		dDeltaCoordsx[i].resize(5 * n[i], 0);
		dDeltaCoordsy[i].resize(5 * n[i], 0);


		for (int j = 0; j < numPeaks[i]; j++){
			for (int k = 0; k < numMasses; k++){
				int previous_x = previousMatch[j*numMasses + k].x;
				int previous_y = previousMatch[j*numMasses + k].y;
				dpreviousMatchx[i][j*numMasses + k] = previous_x;
				dpreviousMatchy[i][j*numMasses + k] = previous_y;


			}
		}

		device_D1[i].resize(n[i], 0);
		device_M1[i].resize(n[i], 0);
	}

	//cout << "All right herhe" << endl;

	vector<vector<int>> diff(total);
	vector<vector<int>> massInt(total);
	vector<vector<int>> peakInt(total);
	vector<vector<int>> deltaInt(total);
	vector<thrust::device_vector<int>> d_diff(total);
	vector < thrust::device_vector<int>> ddeltaInt(total);

	for (int i = 0; i < total; i++){

		int numblocks = 32;
	    int numdeltathreads = n[i] / numblocks + 1;

		diff[i].resize(n[i]);
		convert(cutMasses, massInt[i], cutPeaksMonoMass[i], peakInt[i], ptm, deltaInt[i]);
		calcDiff(massInt[i], peakInt[i], diff[i]);

		d_diff[i]=diff[i];
		ddeltaInt[i]=deltaInt[i];
		//cudaMemcpyAsync(thrust::raw_pointer_cast(&d_diff[i][0]), &diff[i][0], diff[i].size()*sizeof(int), cudaMemcpyHostToDevice, stream[i]);
		//cudaMemcpyAsync(thrust::raw_pointer_cast(&ddeltaInt[i][0]), &deltaInt[i][0], deltaInt[i].size()*sizeof(int), cudaMemcpyHostToDevice, stream[i]);


	/*	cudaEvent_t device_start0, device_stop0;
		cudaEventCreate(&device_start0);
		cudaEventCreate(&device_stop0);
		cudaEventRecord(device_start0, NULL);*/
		KernelInitialDeltaDiag2 << <numblocks, numdeltathreads ,0,stream[i]>> >(numMasses, numPeaks[i], thrust::raw_pointer_cast(&d_diff[i][0]), thrust::raw_pointer_cast(&dNumDeltaCoords[i][0]), thrust::raw_pointer_cast(&dDeltaCoordsx[i][0]), thrust::raw_pointer_cast(&dDeltaCoordsy[i][0]), 20000 * tolerance, thrust::raw_pointer_cast(&ddeltaInt[i][0]), deltaInt[i].size());
		/*cudaEventRecord(device_stop0, NULL);
		cudaEventSynchronize(device_start0);
		cudaEventSynchronize(device_stop0);
		float msecTotal0 = 0.0f;
		cudaEventElapsedTime(&msecTotal0, device_start0, device_stop0);
		out << " GPU计算delta_diag:" << msecTotal0 << "ms";*/
		//KernelInitialDeltaDiag << <numblocks, numdeltathreads >> >(numMasses, numPeaks, thrust::raw_pointer_cast(&d_diff[0]), thrust::raw_pointer_cast(&dNumDeltaCoords[0]), thrust::raw_pointer_cast(&dDeltaCoordsx[0]), thrust::raw_pointer_cast(&dDeltaCoordsy[0]), 2 * tolerance, thrust::raw_pointer_cast(&ddelta[0]), ddelta.size());

		int internalCellRow = 4;
		int internelCellColumn = 4;
		int blocks = 1, threads = numPeaks[i];
		/*cudaEvent_t device_start, device_stop;
		cudaEventCreate(&device_start);
		cudaEventCreate(&device_stop);
		cudaEventRecord(device_start, NULL);*/
		for (int k = 0; k < numSpecific / 2; k++){

			gpu_dp_kernel2 << <blocks, threads,0,stream[i] >> >(numMasses, numPeaks[i], thrust::raw_pointer_cast(&device_M0[i][0]), thrust::raw_pointer_cast(&device_D0[i][0]), thrust::raw_pointer_cast(&device_M1[i][0]), thrust::raw_pointer_cast(&device_D1[i][0]), thrust::raw_pointer_cast(&dpreviousMatchx[i][0]), thrust::raw_pointer_cast(&dpreviousMatchy[i][0]), thrust::raw_pointer_cast(&dNumDeltaCoords[i][0]), thrust::raw_pointer_cast(&dDeltaCoordsx[i][0]), thrust::raw_pointer_cast(&dDeltaCoordsy[i][0]), internalCellRow, internelCellColumn, numMasses / internalCellRow, numMasses%internalCellRow, numPeaks[i] / internelCellColumn);

			gpu_dp_kernel2 << <blocks, threads,0,stream[i] >> >(numMasses, numPeaks[i], thrust::raw_pointer_cast(&device_M1[i][0]), thrust::raw_pointer_cast(&device_D1[i][0]), thrust::raw_pointer_cast(&device_M0[i][0]), thrust::raw_pointer_cast(&device_D0[i][0]), thrust::raw_pointer_cast(&dpreviousMatchx[i][0]), thrust::raw_pointer_cast(&dpreviousMatchy[i][0]), thrust::raw_pointer_cast(&dNumDeltaCoords[i][0]), thrust::raw_pointer_cast(&dDeltaCoordsx[i][0]), thrust::raw_pointer_cast(&dDeltaCoordsy[i][0]), internalCellRow, internelCellColumn, numMasses / internalCellRow, numMasses%internalCellRow, numPeaks[i] / internelCellColumn);

		}
		/*cudaThreadSynchronize();
		cudaEventRecord(device_stop, NULL);
		cudaEventSynchronize(device_start);
		cudaEventSynchronize(device_stop);
		float msecTotal = 0.0f;
		cudaEventElapsedTime(&msecTotal, device_start, device_stop);
		out << " GPU计算6次:" << msecTotal << "ms";
		out << "   no:" << i << ";GPU specific PTM 5: " << device_M1[n - 1] << endl;*/

	}
	clock_t end = clock();
	for (int i = 0; i < total; i++){
		cudaStreamDestroy(stream[i]);
		out << "   no:" << i << ";GPU specific PTM 5: " << device_M1[i][n[i] - 1] << endl;
	}
	out << "--------------------------------------------" << endl;
	out << "--------------------------------------------" << endl;
	out << "Total time:" << end - start << "ms" << endl;
}

bool copare( pair<int,int>& pfirst, pair<int,int>& psecond) 
{                                                    
 return pfirst.second>psecond.second;
}
struct cmp
{
	bool operator()(pair<int, int> &lhs, pair<int, int> &rhs)
	{
		return lhs.second > rhs.second;
	}
};
void filterProtein(vector<vector<mass_t>>& cutMasses,vector<MonoSpectrum>& spectrum,vector<vector<mass_t>>& cutPeaksMonoMass,vector<vector<CutPeak>>& cutPeaks){
	clock_t msfilter_start = clock();
	char *path = "output\\time_result.txt";
	ofstream out(path);
	int limit = 50;
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
	
	vector<vector<pair<int,int>>> diagScore(spectraSize);
	for (int i = 0; i<spectraSize ; i++){
		for (int j = 0; j < index[i].size(); j++){
			int proteinIndex = index[i][j];
			//int value=approxConvolution2(cutMasses[index[i][j]], cutPeaksMonoMass[i], 0.015);
			int value=msFilter(cutMasses[index[i][j]], cutPeaksMonoMass[i], 0.015);
			pair<int, int> tempair;
			tempair = make_pair(proteinIndex, value);
			diagScore[i].push_back(tempair);
		}
		sort(diagScore[i].begin(), diagScore[i].end(),cmp());
	}
	out << "here" << endl;
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

	int nump = 0;
	map<int, int> resu;
	for (int i = 0; i < 60; i++){
		resu.insert(pair<int, int>(i, 0));
	}
	for (int i = 0; i < spectraSize; i++){
		vector<vector<mass_t>> tempMass;
		if (diagScore[i].size()>0){
			int minIndex = 10<diagScore[i].size() ? 10 : diagScore[i].size();
			for (int k = 0; k < minIndex; k++){
				//vector<mass_t>tmp;
				int p = diagScore[i][k].first;
				//for (int j = 0; j < cutMasses[p].size(); j++){
				//	tmp.push_back(cutMasses[p][j]);
				//}
				tempMass.push_back(cutMasses[p]);
				//out << diagScore[i][k].first << "   ";
			}
			//out << endl << "----" << endl;;
		}
		
		vector<int> tscore;
		gpuFunction(tempMass, cutPeaks[i], cutPeaksMonoMass[i], ptm,tscore);
		
		int maxelement = *max_element(tscore.begin(), tscore.end());
		if (maxelement<60)
			resu[maxelement]++;
		if (maxelement > 40){
			out << i << ":" << maxelement << endl;
			if (maxelement > 40)
				nump++;
		}
		
	}

	clock_t dp_end = clock();
	out << "the identified protein num(>30):" << nump<<endl;
	out << "dp time:" << dp_end - msfilter_end << "ms" << endl;
	out<<"dp time:" << (dp_end - msfilter_end)/60000 << "min" << endl;
	out << "--------";
	for (int i = 0; i < 60; i++)
		out << i << " :" << resu[i] << endl;
	
	//cout << "过滤出:"<<index.size() << endl;

}
__global__ void testt(int m){
	int idx = blockDim.x*blockIdx.x + threadIdx.x;
	printf("idx: %d", m);
}
__device__ bool findVal(int* a, int value,int m){
	bool result = false;
	for (int i = 0; i < m; i++){
		if (a[i] == value){
			result = true;
			break;
		}
	}
	return result;
}
__device__ void BubbleSort(double *array,int n)  
{  
	for (int i = 0; i<n - 1; i++)
	{
		for (int j = n - 1; j>i; j--)
		{
			if (array[j] < array[j - 1]){
				int tmp = array[i];
				array[i] = array[j];
				array[j] = tmp;
			}

		}
	}
}  
__global__ void gpuFilter(){
	 int *c=(int *)malloc(sizeof(int)*2000000);
	//int *c = new int[200000];
	 int temp = 0;
	 int t = 5012;
	 for (int i = 1999000; i < 2000000; i++){
		 c[i] = temp++;
		 c[i] = temp*t;
	 }
	 free(c);
	//printf("进来了!");
}
__global__ void gpuFilter(int* cutMass_begin, int* cutMass_end, double* a,double* b,int m,int proteins){
//__global__ void gpuFilter(){
	int idx = blockDim.x*blockIdx.x + threadIdx.x;	
	if (idx < proteins){
		//first convolution stage
		double amino_acid_min = 75.06;
		//clock_t c_start = clock();
		int t = 5;
		int index_start = cutMass_begin[idx];
		int index_end = cutMass_end[idx];
		int n = index_end - index_start + 1;
		double ma = a[index_end];
		double mb = b[m - 1];
		//double ee = 1;xl
		int maxee = (ma + mb) +2;

		//	thrust::device_vector<double> begin;
		//	thrust::device_vector<double> end;
		//	thrust::device_vector<double> diff;
	    double begin[100];
		double end[100];
		double diff[100];
		int size = 0;

		//vector<int> index;
		//	thrust::device_vector<int> tempd;
		//	thrust::device_vector<int> cc(maxee+2, 0);
		int tempd[100];
		 int *cc=(int *)malloc(sizeof(int)*maxee);
		//for (int i = 0; i < 30000; i++)
		//	cc[i] = 0;

		for (int i = index_start; i < index_end + 1; i++){         //calculate the convolution array at resolution ee>e
			double protmass = a[i] + mb;
			for (int j = 0; j < m; j++){
				int d = protmass - b[j];
				cc[d] ++;
				if (cc[d]>t){
					//find function
					if (!findVal(cc, d, size)){
						//index.push_back(j*n + i);
						tempd[size] = d;
						diff[size] = a[i]-b[j];
						size++;
					}
				}
			}
		}
		//cout << "1Da num:" << index.size() << endl;

		//stage2
		int stage2size = 0;
		if (size > 100)
			size = 100;
		BubbleSort(diff, size);
		double e = 0.025;
		for (int i = 0; i < size; i++){
			if (i + 1 < size&&diff[i + 1] - diff[i] + 2  < amino_acid_min){
				begin[stage2size] = (diff[i] - 1);
				end[stage2size] = (diff[i + 1] + 1);
				stage2size++;
				i++;
			}
			else{
				begin[stage2size] = (diff[i] - 1);
				end[stage2size] = (diff[i + 1] + 1);
				stage2size++;
			}
		}

		t = 0;
		//int max = (ma + mb) / e+2;
		//thrust::device_vector<int> c(max+2, 0);
		// int *c=(int *)malloc(sizeof(int)*max);

		for (int k = 0; k < stage2size; k++){

			int i = 0, j = 0;
			double first = begin[k];
			double last = end[k];
			int limit = ((last - first) / e) + 2;
		    int *c = new int[limit];
			//cout << diff[index[k]] << " ";

			while (i < n&&j < m){
				double d = a[index_start + i]-b[j];
				if (d >= last){
					j++;
				}
				else{
					i++;
					if (d > first){
						//  cout << d << endl;
						int indexd = (d -first) / e;
						c[indexd] ++;

					}

				}
			}
			int sum = c[limit - 1] + c[limit - 2];
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

		//printf("CUDA msFilter: %d", t);
	}
	

	
}
__global__ void gpuApproxConvolution(int* cutMass_begin, int* cutMass_end, double* a, double* b, int m, int proteins){
	//int idx = blockDim.x*blockIdx.x + threadIdx.x;
	int idx = 2;
	//if (idx == 3){
		int index_start = cutMass_begin[idx];
		int index_end = cutMass_end[idx];
		int n = index_end - index_start + 1;
		double ma = a[index_end];
		double mb = b[m - 1];

		int score = 0;
		double e = 0.015;
		int max = ((ma+mb)/e)+2;
		//int *c = (int *)malloc(sizeof(int)*max);
		int *c = new int[max];
		//memset(c, 0, sizeof(int)*max);
		//vector<int> c(max + 2, 0);
		//double dmax = 0;
		//testt << <1, 1 >> >(score);
		for (int i = index_start; i <= index_end; i++){
			for (int j = 0; j < m; j++){
				int d = (b[j] - a[i] + ma) / e;
				//int d = 20;
				c[d]++;
				if (c[d]>score){
					score = c[d];
				}
				/*
				d++;
				c[d] = c[d] + 1;
				if (c[d] > score){
					score = c[d];
				}
				*/
			}
		}
		
		free(c);
		//printf("CUDA execution time: %d", score);
	//}
	//__syncthreads();
}
__global__ void calScore(double* a, double* b, int n, int m, double ma, double mb, double e,int* c,int num,int* score,int max){

	int idx = blockDim.x*blockIdx.x + threadIdx.x;
	if (idx < num){
		int i = idx%n;
		int j = idx / n;
		int d = ( a[i]-b[j] + mb)/0.015 ;
		//c[d]++;
		c[d]=atomicAdd(&(c[d]), 1);
		
		if (c[d]>score[0]){
			atomicExch(&(score[0]), c[d]);
			//score = c[d];
			//printf("%d   \n", score[0]);
		}
		

	}
	

}
//有毒
__global__ void gpuExactConvolution(int* cutMass_begin, int* cutMass_end, double* a, double* b, int m, int proteins){
	int idx = blockDim.x*blockIdx.x + threadIdx.x;
	//int idx = 2;
	int index_start = cutMass_begin[idx];
	int index_end = cutMass_end[idx];
	int n = index_end - index_start + 1;
	//vector<double> c(m*n, 0);
	//double *c = (double *)malloc(sizeof(double)*(m*n));
	int num = m*n;
	double *c = new double[num];
	int size = 0;
	for (int i = index_start; i <= index_end; i++){         //calculate the convolution array at resolution ee>e
		for (int j = 0; j < m; j++){
			c[size]=(b[j] - a[i]);
			size++;
		}
	}
	
	//BubbleSort(c, size);
	double e = 0.015;
	//sort(c.begin(), c.end());
	//cout << c[c.size()-1] << endl;
	int score = 0;
	int k = 0;
	for (int i = 0; i < size; i++){
		k = 0;
		double p = c[i];
		if (p != 0){
			for (int j = i + 1; j < size && (j - i) < 50; j++){
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
	free(c);
}

void approxConvolutionmy(const vector<double>& masses, const vector<double>& peaks, double e){
	
	int n = masses.size();
	int m = peaks.size();
	double mb = peaks[m - 1];
	double ma = masses[n - 1];
	int score = 0;
	int max = (ma + mb) / e +2;
	//vector<int> c(max+2, 0);
	clock_t start = clock();
	int *c = new int[max];
	memset(c,0,(max)*sizeof(int));

	int *dev_c = 0;
	cudaMalloc((void**)&dev_c, max * sizeof(int));
	
	for (int i = 0; i < n; i++){
		for (int j = 0; j < m; j++){
			int d = (peaks[j] - masses[i]+ma) / e;
			//cout << "d::" << d <<" m-j:"<<m-j<< endl;
			c[d] = c[d] + 1;
			if (c[d]>score){
				score = c[d];
			//	dmax = peaks[j] - masses[i];
			}
			d++;
			c[d] = c[d] + 1;
			if (c[d]>score){
				score = c[d];
			}
		}
	}
	
	//return score;
	clock_t end = clock();
	std::cout << "一次计算分数时间(approx版本2):" << end - start << "ms" << endl;
	//cout << "p值:" << dmax << endl;
	std::cout << "score:" << score << endl;

}