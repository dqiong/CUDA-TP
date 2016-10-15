
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <cuda.h>
#include <device_functions.h>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/sort.h>
#include <thrust/generate.h>
#include <thrust/copy.h>

#include <stdlib.h>
#include <cstdlib>
#include <time.h>
#include <stdio.h>
#include <string.h>

#include "FileHandle.h"
#include "FilterProtein.cuh"
#include "Handle.cuh"

using namespace std;
#define CUDA_CALL(x){const cudaError_t a=(x);if(a!=cudaSuccess){printf("\n Cuda Error: %s (err_num=%d) \n",cudaGetErrorString(a),a);cudaDeviceReset();assert(0);}}

int main(){
	//gpuFunction2();
	
	int ppm = 15;
	vector<ProteinSequence> proteins;
	Config config;
	char sequence_path[512];
	char spectrum_path[512];
	strcpy(sequence_path, "data\\human.fasta");
	strcpy(spectrum_path, "data\\H4_result.msalign");

	//deal with proteins database
	config.init_with_defaults();
	FileHandle fh;
	fh.readProteins(sequence_path, &config, proteins);
	vector<vector<mass_t>> cutMasses(proteins.size());
	for (int i = 0; i < proteins.size(); i++){
		proteins[i].calcNonModifiedCutMasses(cutMasses[i]);
	}
	//deal with spectrum
	vector<FragmentSet> fragments;
	vector<MonoSpectrum> spectrum;
	vector<vector<CutPeak>> cutPeaks;
	fh.readAllPeaks(spectrum_path, spectrum, fragments, cutPeaks, ppm);
	vector<vector<mass_t>> cutPeaksMonoMass(cutPeaks.size());
	for (int i = 0; i < cutPeaks.size(); i++){
		for (int j = 0; j < cutPeaks[i].size(); j++){
			cutPeaksMonoMass[i].push_back(cutPeaks[i][j].monoMass);
		}
	}
	/* 蛋白质
	char *path = "output\\numdis.txt";
	ofstream out(path);
	int pmax = 0;
	for (int i = 0; i < proteins.size(); i++){
		if (proteins[i].aminoAcidString.length()>pmax)
			pmax = proteins[i].aminoAcidString.length();
	}
	cout << pmax << endl;
		map<int, int> pr;
		int psize = 800;
	for (int i = 0; i < psize+4;i++){
		pr.insert(pair<int, int>(i, 0));
	}
	for (int i = 0; i < proteins.size(); i++){
		int indx = proteins[i].aminoAcidString.length();
		pr[indx]++;
	}
	vector<int> indexp;
	for (int i = 0; i < psize;i+=4){
		indexp.push_back(i + 4);
		 pr[i] = pr[i] + pr[i + 1]+ pr[i + 3] + pr[i + 2];
		out << i+4 << " ";
		
	}
	cout << indexp.size() << endl;
	out << endl;

	for (int i = 0; i < psize;i+=4){
		out << pr[i] << " ";

	}
	*/
	/*谱图
	map<int, int> spect;
	for (int i = 0; i < 200;i++){
		spect.insert(pair<int, int>(i, 0));
	}
	for (int i = 0; i < cutPeaks.size(); i++){
		int indx = cutPeaks[i].size() / 2;
		spect[indx]++;
	}
	for (int i = 0; i < 200;i++){
		out << i << " ";
	}
	out << endl;
	int numbig = 0;
	for (int i = 0; i < 200;i++){
		if (i>20&&i<80)
			numbig += spect[i];
		out << spect[i] << " ";
	}
	cout << numbig<<endl<<"最大比重：" << (double)numbig / 1245 << endl;
	*/
	filterProtein(cutMasses, spectrum, cutPeaksMonoMass,cutPeaks);
	
	
	
	return 0;
}

