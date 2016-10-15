#ifndef __FILEHANDLE_H__
#define __FILEHANDLE_H__

#include "includes.h"
#include "MonoSpectrum.h"
#include "MonoFragments.h"
#include "ProteinSequence.h"
using namespace std;

class FileHandle{

public:
	bool FileHandle::readProteins(char *filePath, Config *_config,vector<ProteinSequence>& proteins)
	{

		ifstream inputStream(filePath);

		if (!inputStream.good())
		{
			cout << "Error: couldn't open file for reading :" << endl << filePath << endl;
			end(1);
		}

		char buff[4096];
		int markNumber = 0;
		/*
		if (inputStream.gcount() < 1)
		{
			cout << "Warning: no sequence was read from :" << endl << filePath << endl;
			return false;
		}*/
		while (inputStream.getline(buff, 4096)){
			ProteinSequence ps;
			ps.config = _config;
			ps.initial();
			// check for fasta name
			if (buff[0] == '>'){
				ps.setFastaName(buff);
				ps.setId(markNumber);
				markNumber++;
				//bool skipReadLine = true;
				while (1){
					//skipReadLine = false;
					//if (!skipReadLine){
						inputStream.getline(buff, 1024);
					//}
					int inputLength = inputStream.gcount();
					Peptide peptide;
					peptide.parse_from_string(ps.config, string(buff));
					//cout << string(buff) << endl;
					const vector<int>& lineAminoAcids = peptide.get_amino_acids();
					unsigned int a;
					for (a = 0; a < lineAminoAcids.size(); a++)
						ps.aminoAcids.push_back(lineAminoAcids[a]);

					ps.aminoAcidString += peptide.as_string(ps.config);

					if (inputStream.eof() || inputStream.gcount() <= 1)
						break;

				}
				proteins.push_back(ps);
			}

		}
		cout << "proteins numbers:" << markNumber << endl;
		return true;
	}

	void readAllPeaks(char *spectrumPath, vector<MonoSpectrum>& spectrum,vector<FragmentSet>& fragments,vector<vector<CutPeak>>& cutPeaks,mass_t ppm){

		

		ifstream ifs(spectrumPath);
		if (!ifs.good())
		{
			cout << "Error: couldn't open mono isotopic mass list for reading:" << endl;
			cout << spectrumPath << endl;
			end(1);
		}


		char buff[128];
		int numPeaks = 0;
		while (ifs.getline(buff, 128))
		{
			char* id_str = strstr(buff, "ID");   //read id
			if (id_str){
				MonoSpectrum mosp;
				FragmentSet fragment;
				vector<CutPeak> cutpeak;
				//mosp.setParentPeakPPM(ppm);
				//mosp.setFragmentPeakPPM(ppm);
				mosp.spectrumIdx = numPeaks;
				numPeaks++;
				mosp.parentMonoMass = -1;
				mosp.totalIntensity = -1;

				if (ifs.getline(buff, 128)){   //read precursor_mass;
					char* precurMass_str = strstr(buff, "PRECURSOR_MASS");
					if (precurMass_str){
						char* equal_token = strstr(buff, "=");
						equal_token++;
						mosp.parentMonoMass = atof(equal_token);
						if (mosp.parentMonoMass <= 0 || mosp.parentMonoMass > 1E7){
							cout << "Error bad value parsed for precursor mass: " << mosp.parentMonoMass << endl;
							cout << "In line: " << buff << endl;
							end(1);
						}
					}
					mass_t parentPPMShift = (ppm * 1E-6 * mosp.parentMonoMass);
					mosp.minParentMonoMass = mosp.parentMonoMass - parentPPMShift;
					mosp.maxParentMonoMass = mosp.parentMonoMass + parentPPMShift;
				}

				score_t totalReadIntensity = 0;
				mosp.monoPeaks.clear();
				mass_t scaleFactor = ppm * 1E-6;
				while (ifs.getline(buff, 64)){
					if (string(buff).length() == 0){
						break;
					}
					MonoPeak peak;
					istringstream iss(buff);
					//iss >> peak.monoMass >> peak.intensity >> peak.orgMz >> peak.charge;
					iss >> peak.monoMass >> peak.intensity;
					if (peak.monoMass <= 0)
						continue;

					mass_t maxPeakShift = peak.monoMass * scaleFactor;

					peak.minMass = peak.monoMass - maxPeakShift;
					peak.maxMass = peak.monoMass + maxPeakShift;

					if (peak.intensity < 0)
						peak.intensity = 1;

					totalReadIntensity += peak.intensity;

					mosp.monoPeaks.push_back(peak);
				}
				sort(mosp.monoPeaks.begin(), mosp.monoPeaks.end());
				// normalize intensity
				if (mosp.totalIntensity <= 0){
					mosp.totalIntensity = totalReadIntensity;
				}
				else{
					score_t normFactor = mosp.totalIntensity / totalReadIntensity;
					unsigned int i;
					for (i = 0; i < mosp.monoPeaks.size(); i++)
						mosp.monoPeaks[i].intensity *= normFactor;
				}

				// enter peaks into map
				mosp.monoPeaksMap.resize(mosp.monoPeaks.size() + 2);
				mosp.monoPeaksMap[0] = -1;
				mosp.monoPeaksMap[mosp.monoPeaks.size() + 1] = 1E9;
				unsigned int i;
				for (i = 0; i < mosp.monoPeaks.size(); i++)
					mosp.monoPeaksMap[i + 1] = mosp.monoPeaks[i].monoMass;

				fragment.initECD();
				spectrum.push_back(mosp);
				fragments.push_back(fragment);
				mosp.convertToCutPeaks(fragment, cutpeak);
				cutPeaks.push_back(cutpeak);

			}
				
		}

		cout << "number of spectra:" << spectrum.size() << endl;
		/*
		for (int i = 0; i < cutPeaks[0].size(); i++){
			cout << cutPeaks[0][i].monoMass << endl; 
		}
		*/

	}



};
#endif