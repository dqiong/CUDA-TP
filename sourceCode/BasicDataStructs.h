
#ifndef __BASICDATASTRUCTS_H__
#define __BASICDATASTRUCTS_H__

#include "includes.h"

class Config;

//多肽类,一个多肽是蛋白质的一个片段,一串字符串
class Peptide {
private:

	mass_t mass; // mass of aas + terminals, doesn't include 19
	int n_term; // N_TERM or a modification (not the aa N-terminal to the peptide)
	int c_term; // C_TERM or a modification (not the aa C-terminal to the peptide)

	vector<int>    amino_acids;
	vector<mass_t> gaps;     
public:
	Peptide() : mass(0), n_term(N_TERM), c_term(C_TERM) {};

	void clear() { mass = 0; n_term = N_TERM; c_term = C_TERM; amino_acids.clear(); gaps.clear(); }

	bool operator== (const Peptide& other) const
	{
		if (amino_acids.size() != other.amino_acids.size())
			return false;
		unsigned int i;
		for (i = 0; i<amino_acids.size(); i++)
		if (amino_acids[i] != other.amino_acids[i])
			return false;
		return true;
	}

	mass_t get_mass() const { return mass; } // mass of aas + terminals, doesn't include 19
	int get_n_term() const { return n_term; }
	int get_c_term() const { return c_term; }
	void set_n_term(int aa) { n_term = aa; }
	void set_c_term(int aa) { c_term = aa; }


	int get_length() const { return amino_acids.size(); } // gaps count as amino acids
	int get_num_aas() const { return amino_acids.size() - gaps.size(); }
	int get_num_gaps() const { return gaps.size(); }

	const vector<int>& get_amino_acids() const { return amino_acids; }

	/*********************************************************************
	Returns the global edit distance between two peptides.
	**********************************************************************/
	float peptide_edit_distance(Config *config, Peptide& other) const;

	// changes the amino acids I->L
	// an Q->K if not at terminal and tolerance > 0.1
	void convert_ILQK(const Config *config);

	// changes all the amino acids to their original form (without PTMs)
	void convert_to_org(const Config *config);

	void convert_IL();

	void reverse();

	string as_string(const Config* config) const;

	void parse_from_string(const Config* config, const string& str);

	void set_peptide(vector<int>& aas, vector<mass_t>& gaps, mass_t mass,
		int n_term_aa = N_TERM, int c_term_aa = C_TERM);

	void set_peptide_aas(const vector<int>& aas) { amino_acids = aas; }

	void calc_mass(const Config *config);

	void calc_expected_breakage_masses(Config *config, vector<mass_t>& break_masses) const;


};



struct CutPeak {
	CutPeak() : monoMass(-1), tolerance(-1), score(0), clusterIdx(-1),
	originalPeakIdx(-1), fragmentInterpertation(-1) {};

	bool operator < (const CutPeak& other) const
	{
		return (monoMass<other.monoMass);
	}

	mass_t	monoMass;
	mass_t	tolerance; // how much +- mass to tolerate with this peak
	score_t score;
	int		clusterIdx;
	int		originalPeakIdx;
	int		fragmentInterpertation; // the idx of the fragment type used to create this mass
};









#endif




