

#include "ProteinSequence.h"
#include "Config.h"

bool ProteinSequence::read(char *filePath, Config *_config, bool verbose)
{
	config = _config;
	ifstream inputStream(filePath);

	if (!inputStream.good())
	{
		cout << "Error: couldn't open file for reading :" << endl << filePath << endl;
		end(1);
	}

	char buff[4096];
	inputStream.getline(buff, 4096);
	if (inputStream.gcount() < 1)
	{
		cout << "Warning: no sequence was read from :" << endl << filePath << endl;
		return false;
	}

	aminoAcids.clear();
	aaModPointers.clear();
	modifications.clear();
	aminoAcidString = "";

	// check for fasta name
	if (buff[0] == '>')
	{
		fastaName = buff + 1;
		inputStream.getline(buff, 4096);
	}

	// check for known mods
	// mods should be of the form a^ a$ aK5 mK25
	if (!strncmp(buff, "#MODS", 4))
	{
		parseModLine(buff, this->modifications);
		inputStream.getline(buff, 4096);
	}

	bool skipReadLine = true;
	while (1)
	{
		if (!skipReadLine)
			inputStream.getline(buff, 4096);
		skipReadLine = false;

		int inputLength = inputStream.gcount();
		Peptide peptide;

		peptide.parse_from_string(config, string(buff));

		const vector<int>& lineAminoAcids = peptide.get_amino_acids();
		unsigned int a;
		for (a = 0; a < lineAminoAcids.size(); a++)
			aminoAcids.push_back(lineAminoAcids[a]);

		aminoAcidString += peptide.as_string(config);

		if (inputStream.eof() || inputStream.gcount() <= 1)
			break;
	}

	if (verbose)
	{
		cout << "Read " << aminoAcids.size() << " :" << endl;
		cout << aminoAcidString << endl;
	}

	//	print();

	return true;
}



mass_t ProteinSequence::calcMass() const
{
	const vector<mass_t>& aa2mass = config->get_aa2mass();
	mass_t totalMass = 0;
	unsigned int i;
	for (i = 0; i < aminoAcids.size(); i++)
		totalMass += aa2mass[aminoAcids[i]];

	return totalMass + totalModMass;
}


void ProteinSequence::calcCutMasses(vector<mass_t>& cutMasses) const
{
	const vector<mass_t>& aa2mass = config->get_aa2mass();
	cutMasses.clear();
	mass_t prefixMass = 0;
	unsigned int i;
	cutMasses.push_back(0);
	for (i = 0; i<aminoAcids.size(); i++)
	{
		prefixMass += aa2mass[aminoAcids[i]];
		if (aaModPointers.size()>0 && aaModPointers[i] >= 0)
		{
			PTM ptm = config->get_PTM(aaModPointers[i]);
			prefixMass += ptm.delta;
		}

		cutMasses.push_back(prefixMass);
	}
}


void ProteinSequence::calcNonModifiedCutMasses(vector<mass_t>& cutMasses) const
{
	const vector<mass_t>& aa2mass = config->get_aa2mass();
	cutMasses.clear();
	mass_t prefixMass = 0;
	unsigned int i;
	cutMasses.push_back(0);
	//cout.setf(ios::showpoint);
	//cout.precision(6);
	//cout.setf(ios::fixed);
	for (i = 0; i < aminoAcids.size(); i++)
	{
		prefixMass += aa2mass[aminoAcids[i]];
		//cout << "  " << aa2mass[aminoAcids[i]] << endl;
		cutMasses.push_back(prefixMass);
	}
}


// reads a line of mods should start with 
void ProteinSequence::parseModLine(char *buff, vector<ModInstance>& mods)
{
	mods.clear();
	istringstream iss(buff);

	while (1)
	{
		string modString = "";
		iss >> modString;

		if (modString.length() < 2)
			break;

		if (modString[0] == '#')
			continue;

		// modstring should be like a,K,12
		// modstring should be like a,^
		unsigned int i;
		for (i = 0; i < modString.length(); i++)
		if (modString[i] == ',' || modString[i] == ':')
			modString[i] = ' ';

		istringstream modStream(modString.c_str());

		string modLabel = "";
		string modAAString = "";
		int modPos = -1;

		modStream >> modLabel >> modAAString >> modPos;

		int modIdx = config->get_PTM_idx(modLabel);
		if (modIdx < 0)
		{
			cout << "Error: unknown mod \'" << modLabel << "\' !" << endl;
			end(1);
		}

		int modAA = config->get_aa_from_label(modAAString);
		if (modAA < 0)
		{
			cout << "Error: unknows mod AA \'" << modAAString << "\' !" << endl;
			end(1);
		}

		ModInstance mod;

		mod.aminoAcid = modAA;
		mod.typeIdx = modIdx;
		mod.pos = modPos;
		modifications.push_back(mod);
	}
}

bool ProteinSequence::checkModificationsConsistency()
{
	aaModPointers.clear();
	aaModPointers.resize(aminoAcids.size(), -1);
	totalModMass = 0;

	size_t i;
	for (i = 0; i < modifications.size(); i++)
	{
		const ModInstance& mod = modifications[i];
		if (mod.aminoAcid == N_TERM)
		{
			if (aaModPointers[0] >= 0)
			{
				cout << "Error: more than one mod on the N-TERM!" << endl;
				end(1);
			}

			aaModPointers[0] = mod.typeIdx;
			totalModMass += config->get_PTM(mod.typeIdx).delta;
			continue;
		}

		if (mod.aminoAcid == C_TERM)
		{
			if (aaModPointers[aminoAcids.size() - 1] >= 0)
			{
				cout << "Error: more than one mod on the N-TERM!" << endl;
				end(1);
			}

			aaModPointers[aminoAcids.size() - 1] = mod.typeIdx;
			totalModMass += config->get_PTM(mod.typeIdx).delta;
			continue;
		}

		if (modifications[i].aminoAcid != aminoAcids[mod.pos])
		{
			cout << "Inconsistency between protein sequence and mod at position " << mod.pos << ":" << endl;
			cout << "mod has " << config->get_aa2label()[mod.aminoAcid] << " sequence has " <<
				config->get_aa2label()[aminoAcids[mod.pos]] << endl;
			end(1);
		}

		// Check that the aa fits the mod
		const PTM& ptm = config->get_PTM(mod.typeIdx);

		unsigned int a;
		for (a = 0; a < ptm.applicableAAs.size(); a++)
		if (aminoAcids[mod.pos] == ptm.applicableAAs[a])
			break;

		if (a == ptm.applicableAAs.size())
		{
			cout << "Error: mod " << ptm.label << " is not applicable to amino acid " <<
				config->get_aa2label()[aminoAcids[mod.pos]] << endl;
			end(1);
		}



		aaModPointers[mod.pos] = mod.typeIdx;
		totalModMass += ptm.delta;

	}



	return true;
}


bool ProteinSequence::applyModLine(char *line)
{
	parseModLine(line, this->modifications);
	checkModificationsConsistency();
	return true;
}

void ProteinSequence::print(ostream& os)
{
	checkModificationsConsistency();
	unsigned int i = 0;
	while (i < aminoAcids.size())
	{
		os << i + 1 << "\t";
		int j = 0;
		while (i < aminoAcids.size() && j < 50)
		{
			os << config->get_aa2label()[aminoAcids[i]];
			i++;
			j++;
		}
		os << "\t" << i << endl;
	}
	os << endl;
	mass_t totalMass = calcMass();
	os << "AA Mass     = " << fixed << setprecision(4) << totalMass - totalModMass << endl;
	os << "Mod Mass    = " << totalModMass << endl;
	os << "Total + H2O = " << totalMass + MASS_H2O << endl;
}


void ProteinSequence::printPrefixCutMasses(int numCuts) const
{
	const vector<string>& aa2label = config->get_aa2label();
	vector<mass_t> cutMasses;

	calcCutMasses(cutMasses);

	unsigned int i;
	for (i = 0; i < cutMasses.size() && i < (unsigned)numCuts; i++)
	{
		cout << i + 1 << "\t";
		cout << setprecision(4) << fixed << cutMasses[i + 1] << "\t";
		unsigned int j;
		for (j = 0; j <= i; j++)
			cout << aa2label[aminoAcids[j]];
		cout << endl;
	}
}

