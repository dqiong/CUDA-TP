


#ifndef __MSMSSPECTRUM_H__
#define __MSMSSPECTRUM_H__

#include "includes.h"
#include "ProteinSequence.h"
#include "BasicDataStructs.h"
#include "MonoFragments.h"



typedef enum AlignmentDir { ALIGN_DIR1, ALIGN_DIR2, ALIGN_BOTH } AlignmentDir;

class Spectrum;


static int runningSpectrumIdx = 0;



// A peak that is read from the spectrum file
struct MonoPeak {
	MonoPeak() : monoMass(0), minMass(0), maxMass(0), orgMz(0), charge(0), intensity(0) {};

	bool operator< (const MonoPeak& other) const
	{
		return (monoMass < other.monoMass);
	}

	mass_t monoMass;
	mass_t minMass; // lower boundry according to fragmentPeakPPM
	mass_t maxMass; // upper boundry according to fragmentPeakPPM

	mass_t orgMz;
	int    charge;

	intensity_t intensity;

};





// A spectrum that is read from the spectrum files
// Same class for all MS levels
class MonoSpectrum {
	friend class SpectrumForest;
	friend class FileHandle;
public:
	MonoSpectrum() : spectrumPath(NULL), depth(-1), spectrumIdx(-1), parentSpectrumIdx(-1),
		parentPeakIdx(-1), parentPeakPPM(-1), fragmentPeakPPM(-1),
		parentMonoMass(-1), minParentMonoMass(-1), maxParentMonoMass(-1),
		totalIntensity(0) {};

	char	*getSpectrumPath() const { return spectrumPath; }

	int		getDepth() const { return depth; }
	void	setDepth(int _depth) { depth = _depth; }

	int		getSpectrumIdx() const { return spectrumIdx; }
	int		getParentSpectrumIdx() const { return parentSpectrumIdx; }

	mass_t  getParentMonoMass() const { return parentMonoMass; }
	void	setParentMonoMass(mass_t pm) { parentMonoMass = pm; }
	mass_t	getParentPeakPPM() const { return parentPeakPPM; }
	void	setParentPeakPPM(mass_t _parentPPM) { parentPeakPPM = _parentPPM; }

	mass_t  getFragmentPeakPPM() const { return fragmentPeakPPM; }
	void	setFragmentPeakPPM(mass_t _fragmentPPM) { fragmentPeakPPM = _fragmentPPM; }

	score_t getTotalIntensity() const { return totalIntensity; }

	int		getNumMonoPeaks() const { return monoPeaks.size(); }
	const vector<MonoPeak>& getMonoPeaks() const { return monoPeaks; }

	int		getNumDaughterSpectra() const { return daughterSpectraIdxs.size(); }
	const vector<int>& getDaughterSpectra() const { return daughterSpectraIdxs; }

	int		readMonoMassPeakList(char *_spectrumPath, mass_t _precursorPeakPPM,
		mass_t _fragmentPeakPPM = -1);

	int		readAllPeakList(char *_spectrumPath, mass_t _precursorPeakPPM,
		mass_t _fragmentPeakPPM = -1);

	int		getNumMonoPeaksInRange(mass_t minMass, mass_t maxMass, int *firstIdx) const;

	void	findFragmentPeaks(const ProteinSequence& ps, const FragmentSet& fs,
		vector< vector<int> >& peakIdxs) const;

	intensity_t	scoreSpectrumMatch(const ProteinSequence& ps, const FragmentSet& fs,
		int& numMatchedPeaks) const;

	void	convertToCutPeaks(const FragmentSet& fs, vector<CutPeak>& cutPeaks) const;

	void	printMatchedFragPositions(const FragmentSet& fs, const ProteinSequence& ps,
		vector< vector<int> >& peakIdxs) const;

	void	print(ostream& os = cout) const;

	//void setParentPeakPPM(mass_t ppm){ parentPeakPPM = ppm; }
	//void setFragmentPeakPPM(mass_t ppm){ fragmentPeakPPM = ppm; }

private:

	char * spectrumPath;

	int	   depth;				// depth in the tree MS1/MS2/MS3...  ,no use
	int	   spectrumIdx;
	int	   parentSpectrumIdx;
	int	   parentPeakIdx;

	mass_t parentPeakPPM;
	mass_t fragmentPeakPPM;

	mass_t parentMonoMass;
	mass_t minParentMonoMass;
	mass_t maxParentMonoMass;

	score_t totalIntensity;

	vector<MonoPeak> monoPeaks;
	vector<int>      daughterSpectraIdxs;

	vector<mass_t> monoPeaksMap;
};





#endif

