
#include "MonoSpectrum.h"


/*****************************************************************
Reads the monisotopic mass from a peak list.
- Assumes that the mono isotopic parent mass (and intensity - otional)
appear in the first line that starts with numbers otherwise.
- Assumes that each line is a peak mass and intensity.
if the intensity is missing, each peak is assigned an intensity 1.

returns the number of peaks read (-1 if there is an error).
******************************************************************/
int		MonoSpectrum::readMonoMassPeakList(char *_spectrumPath,
	mass_t _precursorPeakPPM, mass_t _fragmentPeakPPM)
{
	parentPeakPPM = _precursorPeakPPM;
	fragmentPeakPPM = _fragmentPeakPPM;
	if (_fragmentPeakPPM<0)
		fragmentPeakPPM = parentPeakPPM;


	spectrumPath = _spectrumPath;
	ifstream ifs(spectrumPath);
	if (!ifs.good())
	{
		cout << "Error: couldn't open mono isotopic mass list for reading:" << endl;
		cout << spectrumPath << endl;
		end(1);
	}


	char buff[128];
	while (ifs.getline(buff, 128))
	{
		if (buff[0] == '#')
			continue;
		break;
	}

	parentMonoMass = -1;
	totalIntensity = -1;
	istringstream iss(buff);
	iss >> parentMonoMass >> totalIntensity;
	if (parentMonoMass <= 0 || parentMonoMass>1E7)
	{
		cout << "Error bad value parsed for precursor mass: " << parentMonoMass << endl;
		cout << "In line: " << buff << endl;
		end(1);
	}

	mass_t parentPPMShift = (parentPeakPPM * 1E-6 * parentMonoMass);
	minParentMonoMass = parentMonoMass - parentPPMShift;
	maxParentMonoMass = parentMonoMass + parentPPMShift;

	score_t totalReadIntensity = 0;
	monoPeaks.clear();

	mass_t scaleFactor = fragmentPeakPPM * 1E-6;
	while (ifs.getline(buff, 64))
	{
		MonoPeak peak;

		istringstream iss(buff);
		iss >> peak.monoMass >> peak.intensity;

		if (peak.monoMass <= 0)
			continue;

		mass_t maxPeakShift = peak.monoMass * scaleFactor;

		peak.minMass = peak.monoMass - maxPeakShift;
		peak.maxMass = peak.monoMass + maxPeakShift;

		if (peak.intensity<0)
			peak.intensity = 1;

		totalReadIntensity += peak.intensity;

		monoPeaks.push_back(peak);
	}

	// normalize intensity
	if (totalIntensity <= 0)
	{
		totalIntensity = totalReadIntensity;
	}
	else
	{
		score_t normFactor = totalIntensity / totalReadIntensity;
		unsigned int i;
		for (i = 0; i<monoPeaks.size(); i++)
			monoPeaks[i].intensity *= normFactor;
	}


	// enter peaks into map
	monoPeaksMap.resize(monoPeaks.size() + 2);
	monoPeaksMap[0] = -1;
	monoPeaksMap[monoPeaks.size() + 1] = 1E9;
	unsigned int i;
	for (i = 0; i<monoPeaks.size(); i++)
		monoPeaksMap[i + 1] = monoPeaks[i].monoMass;


	return monoPeaks.size();
}


int		MonoSpectrum::readAllPeakList(char *_spectrumPath, mass_t _precursorPeakPPM,
	mass_t _fragmentPeakPPM)
{
	parentPeakPPM = _precursorPeakPPM;
	fragmentPeakPPM = _fragmentPeakPPM;
	if (_fragmentPeakPPM<0)
		fragmentPeakPPM = parentPeakPPM;


	spectrumPath = _spectrumPath;
	ifstream ifs(spectrumPath);
	if (!ifs.good())
	{
		cout << "Error: couldn't open mono isotopic mass list for reading:" << endl;
		cout << spectrumPath << endl;
		end(1);
	}


	char buff[128];
	while (ifs.getline(buff, 128))
	{
		if (buff[0] == '#')
			continue;
		break;
	}

	parentMonoMass = -1;
	totalIntensity = -1;
	istringstream iss(buff);
	iss >> parentMonoMass >> totalIntensity;
	if (parentMonoMass <= 0 || parentMonoMass>1E7)
	{
		cout << "Error bad value parsed for precursor mass: " << parentMonoMass << endl;
		cout << "In line: " << buff << endl;
		end(1);
	}

	mass_t parentPPMShift = (parentPeakPPM * 1E-6 * parentMonoMass);
	minParentMonoMass = parentMonoMass - parentPPMShift;
	maxParentMonoMass = parentMonoMass + parentPPMShift;

	score_t totalReadIntensity = 0;
	monoPeaks.clear();

	mass_t scaleFactor = fragmentPeakPPM * 1E-6;
	while (ifs.getline(buff, 64))
	{
		MonoPeak peak;

		istringstream iss(buff);
		iss >> peak.monoMass >> peak.intensity >> peak.orgMz >> peak.charge;

		if (peak.monoMass <= 0)
			continue;

		mass_t maxPeakShift = peak.monoMass * scaleFactor;

		peak.minMass = peak.monoMass - maxPeakShift;
		peak.maxMass = peak.monoMass + maxPeakShift;

		if (peak.intensity<0)
			peak.intensity = 1;

		totalReadIntensity += peak.intensity;

		monoPeaks.push_back(peak);
	}

	sort(monoPeaks.begin(), monoPeaks.end());

	// normalize intensity
	if (totalIntensity <= 0)
	{
		totalIntensity = totalReadIntensity;
	}
	else
	{
		score_t normFactor = totalIntensity / totalReadIntensity;
		unsigned int i;
		for (i = 0; i<monoPeaks.size(); i++)
			monoPeaks[i].intensity *= normFactor;
	}


	// enter peaks into map
	monoPeaksMap.resize(monoPeaks.size() + 2);
	monoPeaksMap[0] = -1;
	monoPeaksMap[monoPeaks.size() + 1] = 1E9;
	unsigned int i;
	for (i = 0; i<monoPeaks.size(); i++)
		monoPeaksMap[i + 1] = monoPeaks[i].monoMass;


	return monoPeaks.size();
}



/*****************************************************************************
Takes each mono peak and interperts it as one of the given fragments (c,z, etc)
creates a cut peak for each such iterpertation.
******************************************************************************/
void MonoSpectrum::convertToCutPeaks(const FragmentSet& fs, vector<CutPeak>& cutPeaks) const
{
	const vector<int>& fragmentIdxs = fs.getStrongFragmentIdxs();
	const mass_t toleranceFactor = fragmentPeakPPM * 1E-6;
	const mass_t suffixTolerance = parentMonoMass * parentPeakPPM * 1E-6;

	cutPeaks.clear();
	CutPeak cp;
	cp.monoMass = 0;
	cp.tolerance = suffixTolerance;
	cp.fragmentInterpertation = -1;
	cp.originalPeakIdx = -1;
	cp.score = 0;
	cutPeaks.push_back(cp);

	cp.monoMass = parentMonoMass - MASS_H2O;
	cutPeaks.push_back(cp);

	const double maxCutPeakMass = cp.monoMass;
	unsigned int f;
	for (f = 0; f<fragmentIdxs.size(); f++)
	{
		const int fragIdx = fragmentIdxs[f];
		const MonoFragment& frag = fs.getFragment(fragIdx);
		unsigned int peakIdx;
		for (peakIdx = 0; peakIdx<monoPeaks.size(); peakIdx++)
		{
			const MonoPeak& peak = monoPeaks[peakIdx];
			CutPeak cutPeak;

			const mass_t prefixMass = frag.calcExpCutMass(peak.monoMass, parentMonoMass);
			const mass_t tolerance = (frag.getDirection() == PREFIX ?
				peak.monoMass * toleranceFactor : suffixTolerance);

			if (prefixMass <= 0.0 || prefixMass >= maxCutPeakMass)
				continue;

			cutPeak.monoMass = prefixMass;
			cutPeak.tolerance = tolerance;
			cutPeak.score = 1;
			cutPeak.clusterIdx = 0;
			cutPeak.fragmentInterpertation = fragIdx;
			cutPeak.originalPeakIdx = peakIdx;

			cutPeaks.push_back(cutPeak);
		}
	}

	sort(cutPeaks.begin(), cutPeaks.end());
}


int	MonoSpectrum::getNumMonoPeaksInRange(mass_t minMass, mass_t maxMass, int *firstIdx) const
{
	*firstIdx = -1;

	if (minMass<0)
		minMass = 0;
	if (maxMass>5E8)
		maxMass = 5E8;

	vector<mass_t>::const_iterator location = lower_bound(monoPeaksMap.begin(), monoPeaksMap.end(), minMass);
	while (*location<maxMass)
		++location;
	while (*location>maxMass)
		--location;

	if (*location<minMass)
		return 0;

	int topIdx = location - monoPeaksMap.begin();

	while (*location >= minMass)
		--location;

	++location;
	int botIdx = location - monoPeaksMap.begin();

	*firstIdx = botIdx - 1; // there is a dummy -1 in the begining

	return (topIdx - botIdx + 1);
}



void MonoSpectrum::findFragmentPeaks(const ProteinSequence& ps,
	const FragmentSet& fs, vector< vector<int> >& peakIdxs) const
{
	vector<mass_t> cutMasses;
	ps.calcCutMasses(cutMasses);
	peakIdxs.resize(cutMasses.size());

	const int numFrags = fs.getNumFragments();
	const mass_t ppmMul = fragmentPeakPPM * 1E-6;

	unsigned int cutIdx;
	for (cutIdx = 0; cutIdx<cutMasses.size(); cutIdx++)
	{
		vector<mass_t> fragMasses;
		fs.calcFragmentSetMonoMasses(cutMasses[cutIdx], parentMonoMass, fragMasses);

		peakIdxs[cutIdx].clear();
		peakIdxs[cutIdx].resize(numFrags, -1);
		int f;
		for (f = 0; f<numFrags; f++)
		{
			const mass_t& fragMass = fragMasses[f];
			const mass_t ppmError = ppmMul * fragMass;
			const mass_t minFragMass = fragMass - ppmError;
			const mass_t maxFragMass = fragMass + ppmError;

			int firstPeakIdx = -1;
			int numPeaksInRange = getNumMonoPeaksInRange(minFragMass, maxFragMass, &firstPeakIdx);
			peakIdxs[cutIdx][f] = firstPeakIdx;

			if (numPeaksInRange <= 1)
				continue;

			// choose peak with minimal distance from the expected fragment mass
			int closestIdx = firstPeakIdx;
			mass_t closestPeakDistance = fabs(fragMass - monoPeaks[firstPeakIdx].monoMass);
			int i;
			for (i = 1; i<numPeaksInRange; i++)
			{
				const mass_t peakDistance = fabs(fragMass - monoPeaks[firstPeakIdx + i].monoMass);
				if (peakDistance<closestPeakDistance)
				{
					closestPeakDistance = peakDistance;
					closestIdx = firstPeakIdx + i;
				}
			}
			peakIdxs[cutIdx][f] = closestIdx;
		}
	}
}


/***************************************************************************
Matches a spectrum with a protein sequence.
Sums up the intensity of the matched peaks.
****************************************************************************/
intensity_t	MonoSpectrum::scoreSpectrumMatch(const ProteinSequence& ps,
	const FragmentSet& fs, int& numMatchedPeaks) const
{
	intensity_t totalMatchedIntensity = 0;
	vector< vector<int> > fragPeakIdxs;
	findFragmentPeaks(ps, fs, fragPeakIdxs);

	numMatchedPeaks = 0;

	unsigned int cutIdx;
	for (cutIdx = 0; cutIdx<fragPeakIdxs.size(); cutIdx++)
	{
		unsigned int f;
		for (f = 0; f<fragPeakIdxs[cutIdx].size(); f++)
		{
			if (fragPeakIdxs[cutIdx][f]<0)
				continue;

			numMatchedPeaks++;
			totalMatchedIntensity += monoPeaks[fragPeakIdxs[cutIdx][f]].intensity;
		}
	}
	return totalMatchedIntensity;
}


void MonoSpectrum::print(ostream& os) const
{
	unsigned int i;
	os << fixed << setprecision(3) << this->parentMonoMass << "\t" << this->totalIntensity;
	os << "\t(" << this->minParentMonoMass << "\t" << this->maxParentMonoMass << ")" << endl;
	for (i = 0; i<monoPeaks.size(); i++)
	{
		os << i << "\t" << monoPeaks[i].monoMass << "\t" << monoPeaks[i].intensity << "\t" <<
			"(" << monoPeaks[i].minMass << "\t" << monoPeaks[i].maxMass << ")" << endl;
	}
}


void MonoSpectrum::printMatchedFragPositions(const FragmentSet& fs, const ProteinSequence& ps,
	vector< vector<int> >& peakIdxs) const
{
	vector<mass_t> cutMasses;
	ps.calcCutMasses(cutMasses);

	unsigned int cutIdx;
	for (cutIdx = 0; cutIdx<peakIdxs.size(); cutIdx++)
	{
		unsigned int f;
		for (f = 0; f<peakIdxs[cutIdx].size(); f++)
		{
			if (peakIdxs[cutIdx][f]<0)
				continue;

			const int peakIdx = peakIdxs[cutIdx][f];
			const MonoFragment& frag = fs.getFragment(f);
			const mass_t expPeakMass = frag.calcExpFragMass(cutMasses[cutIdx], parentMonoMass);

			mass_t ppmError = 1E6* ((expPeakMass - monoPeaks[peakIdx].monoMass) / expPeakMass);

			cout << cutIdx << " " << frag.getLabel() << "\t";

			cout << setprecision(3) << fixed << expPeakMass << " " << monoPeaks[peakIdx].monoMass << "\t"
				<< ppmError << endl;
		}
	}
}
