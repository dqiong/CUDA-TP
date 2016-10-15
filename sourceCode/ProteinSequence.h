

#ifndef __PROTEINSEQUENCE_H__
#define __PROTEINSEQUENCE_H__

#include "Config.h"
#include "auxfun.h"
#include "includes.h"

struct ModInstance {
	int  typeIdx;
	int  aminoAcid;
	int  pos;
};

class ProteinSequence {
	friend class FileHandle;
public:
	ProteinSequence() : totalModMass(0), config(NULL), aminoAcidString(""), fastaName("") {};

	void initial(){
		aminoAcids.clear();
		aaModPointers.clear();
		modifications.clear();
		aminoAcidString = "";
	}

	bool read(char *filePath, Config *config, bool verbose = false);

	bool applyModLine(char *line);

	void print(ostream& os = cout);

	mass_t calcMass() const;

	void calcCutMasses(vector<mass_t>& cutMasses) const;

	void calcNonModifiedCutMasses(vector<mass_t>& cutMasses) const;

	const vector<int>& getAminoAcids() const { return aminoAcids; }

	int	  getNumCuts() const { return aminoAcids.size() + 1; }

	void printPrefixCutMasses(int numCuts = 10) const;

	Config * getConfig() const { return config; }

	ProteinSequence& operator = (const ProteinSequence& other)
	{
		aminoAcids = other.aminoAcids;
		aaModPointers = other.aaModPointers;
		totalModMass = other.totalModMass;
		config = other.config;
		aminoAcidString = other.aminoAcidString;
		fastaName = other.fastaName;
		return *this;
	}

	void permuteSequence()
	{
		permuteVector(aminoAcids);
	}
	string getAmnioString(){ return aminoAcidString; }
	string getFastaName(){ return fastaName; }
	void setFastaName(char* c){ fastaName = c; }
	void setId(int id){ proteinId = id; }
	void aminoAcid_push(int i){ aminoAcids.push_back(i); }

	string				aminoAcidString;
private:
	vector<int>         aminoAcids;      // holds only amino acids, no terminal aas
	vector<int>			aaModPointers;

	mass_t				totalModMass;
	Config *			config;

	//this is mark number for protein
	int                 proteinId;

	vector<ModInstance> modifications;

	string				fastaName;


	void parseModLine(char *buff, vector<ModInstance>& mods);

	bool checkModificationsConsistency();

};



#endif
