
#ifndef __CONVERSIONTABLES_H__
#define __CONVERSIONTABLES_H__

#include "includes.h"

class ConversionTables {
public:
	void init_for_standard_aas();

	void add_optional_PTM_aa(int aa, string label, double delta)
	{
		org_aa.push_back(aa);
		aa2mass.push_back(aa2mass[aa] + delta);
		aa2label.push_back(label);
	}

	void make_fixed_mod(int aa, double delta)
	{
		aa2mass[aa] += delta;
		if (aa > 0 && aa <= Val)
			char2mass[aa2char[aa]] += delta;
	}




	// access to 
	int	   get_char2aa(char c) const { return char2aa[c]; }
	mass_t get_char2mass(char c) const { return char2mass[c]; }
	mass_t get_aa2mass(int a) const { return aa2mass[a]; }
	char   get_aa2char(int a) const { return aa2char[a]; }
	string get_aa2label(int a) const { return aa2label[a]; }
	int    get_org_aa(int a) const { return org_aa[a]; }

	const vector<int>& get_char2aa() const { return char2aa; }
	const vector<mass_t>& get_char2mass() const { return char2mass; }
	const vector<mass_t>& get_aa2mass() const { return aa2mass; }
	const vector<char>& get_aa2char() const { return aa2char; }
	const vector<string>& get_aa2label() const { return aa2label; }
	const vector<int>& get_org_aa() const { return org_aa; }

private:
	vector<int>       char2aa;
	vector<mass_t> char2mass;
	vector<mass_t> aa2mass;
	vector<char>      aa2char;
	vector<string>    aa2label;
	vector<int>       org_aa;

};


#endif

