
#include "BasicDataStructs.h"
#include "Config.h"
#include "auxfun.h"


/**********************************************************
Parses a string into amino acids.
Terminals are stored in n_term and c_term.
Final mass = sum mass aas + mass terminals
***********************************************************/
void Peptide::parse_from_string(const Config* config, const string &str)
{
	int i;
	const vector<int>& char2aa = config->get_char2aa();

	n_term = N_TERM;
	c_term = C_TERM;
	amino_acids.clear();
	gaps.clear();

	int str_length = str.length();
	while (str_length > 0 && (str[str_length - 1] == '\n' || str[str_length - 1] == '\t'
		|| str[str_length - 1] == ' '))
		str_length--;

	for (i = 0; i < str_length; i++)
	{
		if (char2aa[str[i]] < 0 && str[i] != '[')
		{
			cout << "Error parsing peptide string!" << endl;
			cout << str << endl << setw(i + 1) << right << "^" << endl;
			cout << "Bad token: " << str[i] << endl;
			end(1);
		}

		// this can only be a single char unmodifed AA
		if (i == str_length - 1)
		{
			int aa_idx = char2aa[str[i]];
			if (aa_idx<c_term || aa_idx>Val)
			{
				cout << "Error parsing peptide string!" << endl;
				cout << str << endl << setw(i + 1) << right << "^" << endl;
				cout << "Bad token: " << str[i] << endl;
				end(1);
			}

			// do not push unmodified terminals
			if (aa_idx > c_term)
				amino_acids.push_back(aa_idx);

			continue;
		}

		// parse gap
		if (str[i] == '[')
		{
			int close_idx = i + 1;
			while (close_idx < str_length && str[close_idx] != ']')
				close_idx++;

			if (close_idx == str_length)
			{
				cout << "Error parsing peptide string!" << endl;
				cout << str << endl << setw(i + 1) << right << "^" << endl;
				cout << "(No terminating ']' found)." << endl;
				end(1);
			}

			double gap_size = NEG_INF;
			istringstream is(str.substr(i + 1, close_idx - i - 1));
			is >> gap_size;

			if (gap_size <= 0)
			{
				cout << "Error parsing peptide string!" << endl;
				cout << str << endl << setw(i + 1) << right << "^" << endl;
				cout << "Bad gap string:" << is.str() << endl;
				end(1);
			}

			gaps.push_back(gap_size);
			amino_acids.push_back(Gap);
			i = close_idx;
			continue;
		}

		// parse an aa label
		// find label end
		{
			int end_idx = i + 1;
			while (end_idx < str_length && char2aa[str[end_idx]]<0 && str[end_idx] != '[')
				end_idx++;

			// look for label
			int label_length = end_idx - i;
			while (label_length>1 && (str[i + label_length - 1] == ' ' || str[i + label_length - 1] == '\t'))
				label_length--;

			string label = str.substr(i, label_length);
			int aa_idx = config->get_aa_from_label(label);
			if (aa_idx < 0)
			{
				cout << "Error parsing peptide string!" << endl;
				cout << str << endl << setw(i + 1) << right << "^" << endl;
				cout << "Bad AA label: \'" << label << "\'" << endl;
				end(1);
			}

			if (label[0] == '^')
			{
				n_term = aa_idx;
				if (amino_acids.size() > 0)
				{
					cout << "Error: N-terminal should start peptide string: " << label << endl;
					end(1);
				}
				i = end_idx - 1;
				continue;
			}
			else if (label[0] == '$')
			{
				c_term = aa_idx;
				break;
			}
			else
			{
				amino_acids.push_back(aa_idx);
				i = end_idx - 1;
			}
			continue;
		}
	}

	//calculate total mass

	int g_count = 0;
	const ConversionTables& session_tables = config->get_session_tables();

	// add N-terminal
	mass = session_tables.get_aa2mass(n_term);

	if (amino_acids.size() > 0)
	{
		// add first aa
		if (amino_acids[0] == Gap)
		{
			mass += gaps[g_count++];
		}
		else
			mass += session_tables.get_aa2mass(amino_acids[0]);

		// add middle aas
		unsigned int i;
		for (i = 1; i < amino_acids.size() - 1; i++)
		{
			if (amino_acids[i] == Gap)
			{
				mass += gaps[g_count++];
			}
			else
				mass += session_tables.get_aa2mass(amino_acids[i]);
		}

		// add last aa
		if (amino_acids.size() > 1)
		{
			int last_idx = amino_acids.size() - 1;
			if (amino_acids[last_idx] == Gap)
			{
				mass += gaps[g_count++];
			}
			else
				mass += session_tables.get_aa2mass(amino_acids[last_idx]);
		}
	}

	// add C-terminal
	mass += session_tables.get_aa2mass(c_term);

}





void Peptide::set_peptide(vector<int>& aas, vector<mass_t>& aa_gaps, mass_t peptide_mass,
	int n_term_aa, int c_term_aa)
{
	amino_acids = aas;
	gaps = aa_gaps;
	mass = peptide_mass;
	n_term = n_term_aa;
	c_term = c_term_aa;
}

// the mass without 19
void Peptide::calc_mass(const Config *config)
{
	const vector<mass_t>& aa2mass = config->get_aa2mass();
	unsigned int i;

	mass = aa2mass[n_term] + aa2mass[c_term];

	for (i = 0; i < amino_acids.size(); i++)
		mass += aa2mass[amino_acids[i]];

	for (i = 0; i < gaps.size(); i++)
		mass += gaps[i];
}



// changes the amino acids I->L
// an Q->K if not at terminal and tolerance > 0.1
void Peptide::convert_ILQK(const Config *config)
{
	unsigned int i;
	for (i = 0; i<amino_acids.size(); i++)
	if (amino_acids[i] == Ile)
		amino_acids[i] = Leu;

	if (config->get_tolerance()>0.1)
	for (i = 0; i < amino_acids.size() - 1; i++)
	if (amino_acids[i] == Lys)
		amino_acids[i] = Gln;
}

// changes the amino acids I->L
void Peptide::convert_IL()
{
	unsigned int i;
	for (i = 0; i < amino_acids.size(); i++)
	if (amino_acids[i] == Ile)
		amino_acids[i] = Leu;
}

void Peptide::reverse()
{
	vector<int> tmp_aas;
	vector<mass_t> tmp_gaps;
	int i;

	tmp_aas.clear();
	for (i = amino_acids.size() - 1; i >= 0; i--)
		tmp_aas.push_back(amino_acids[i]);
	amino_acids = tmp_aas;

	tmp_gaps.clear();
	for (i = gaps.size() - 1; i >= 0; i--)
		tmp_gaps.push_back(gaps[i]);

	gaps = tmp_gaps;
}



string Peptide::as_string(const Config* config) const
{
	const vector<string>& aa2label = config->get_aa2label();
	int g_count = 0;
	unsigned int i;

	string peptide_str = "";

	if (n_term > N_TERM)
		peptide_str += aa2label[n_term];

	for (i = 0; i < amino_acids.size(); i++)
	{
		if (amino_acids[i] == Gap)
		{
			ostringstream os;
			os << gaps[g_count++];

			peptide_str += '[' + os.str() + ']';
		}
		else
			peptide_str += aa2label[amino_acids[i]];
	}

	if (c_term > C_TERM)
		peptide_str += aa2label[c_term];

	return peptide_str;
}



void Peptide::calc_expected_breakage_masses(Config *config, vector<mass_t>& break_masses) const
{
	unsigned int i;
	unsigned int aa_pos = 0;
	int g_count = 0;
	mass_t m = 0;
	break_masses.clear();
	if (amino_acids.size() == 0)
		return;

	const ConversionTables& session_tables = config->get_session_tables();

	break_masses.push_back(0);
	m += session_tables.get_aa2mass(n_term);

	// use pre tables for first amino acid
	if (amino_acids[0] > Gap)
	{
		m += session_tables.get_aa2mass(amino_acids[0]);
		break_masses.push_back(m);
		aa_pos++;
	}

	// use mid tables for center amino acids
	for (i = aa_pos; i < amino_acids.size() - 1; i++)
	{
		if (amino_acids[i] == Gap)
		{
			m += this->gaps[g_count++];
			continue;
		}

		m += session_tables.get_aa2mass(amino_acids[i]);
		break_masses.push_back(m);
	}

	// use suf tables for last amino acids
	if (amino_acids.size() > 1)
	{
		const int aa_pos = amino_acids.size() - 1;
		if (amino_acids[aa_pos] == Gap)
		{
			m += this->gaps[g_count++];
		}
		else
			m += session_tables.get_aa2mass(amino_acids[aa_pos]);

		// add C-terminal
		m += session_tables.get_aa2mass(c_term);
	}

	break_masses.push_back(m);
}



/**************************************************************************
Returns the global edit distance between two peptides.
Gap cost = 1.
d(I,L)=0
d(K,Q)=0
d(F/M*)=0
d(N,D)=0.5;
d(X,X)=0
d(X,Y)=1
***************************************************************************/
float Peptide::peptide_edit_distance(Config *config, Peptide& other_pep) const
{
	int i;
	vector<int> other_aa = other_pep.get_amino_acids();
	const int *pep1 = &amino_acids[0];
	const int pep_len1 = amino_acids.size();
	const int *pep2 = &other_aa[0];
	const int pep_len2 = other_aa.size();
	const int max_width = 5;
	const int ox_met_aa = config->get_aa_from_label(string("M+16"));

	float row1[max_width * 2 + 1], row2[max_width * 2 + 1];
	float *old_row, *new_row;

	if (abs(pep_len1 - pep_len2) >= max_width)
		return static_cast<float>(pep_len1);

	// check that no little switch of two aa can make them the same
	if (pep_len1 == pep_len2)
	{
		int i;
		int err_pos = -1;
		int errs = 0;
		for (i = 0; i<pep_len1; i++)
		{
			if ((pep1[i] != pep2[i]) &&
				!(((pep1[i] == Ile || pep1[i] == Leu) && (pep2[i] == Ile || pep2[i] == Leu)) ||
				((pep1[i] == Gln || pep1[i] == Lys) && (pep2[i] == Gln || pep2[i] == Lys))))
			{
				err_pos = i;
				errs++;
			}
			if (errs>2)
				break;
		}

		if (errs == 0)
			return 0;

		if (errs == 2)
		{
			if ((pep1[err_pos] == pep2[err_pos - 1]) && (pep1[err_pos - 1] == pep2[err_pos]))
				return 1;
		}
	}


	old_row = row1;
	new_row = row2;
	for (i = 0; i < max_width; i++)
	{
		old_row[i] = 9999.0;
		old_row[i + max_width] = static_cast<float>(i);
		new_row[i] = 9999.0;
		new_row[i + max_width] = 9999.0;
	}
	old_row[2 * max_width] = 9999.0;
	new_row[2 * max_width] = 9999.0;

	int start_new_row = max_width - 1;
	for (i = 1; i <= pep_len1; i++)
	{
		int j;

		if (start_new_row > 0)
		{
			new_row[start_new_row] = static_cast<float>(i);
			new_row[--start_new_row] = 9999.0;
		}
		else
			new_row[0] = 9999.0;

		for (j = 1; j <= 2 * max_width; j++)
		{
			int p2_pos = i + j - max_width;
			if (p2_pos<1)
				continue;
			if (p2_pos>pep_len2)
				break;

			float v1, v2, v3, dxy = 0;
			int p1 = pep1[i - 1], p2 = pep2[p2_pos - 1];

			if (p1 != p2)
			{
				dxy = 1;
				if (((p1 == Ile || p1 == Leu) && (p2 == Ile || p2 == Leu)) ||
					((p1 == Gln || p1 == Lys) && (p2 == Gln || p2 == Lys)))
				{
					dxy = 0;
				}
				else if (((p1 == Asn && p2 == Asp) || (p1 == Asp && p2 == Asn)) ||
					((pep1[i] == Gln && p2 == Glu) || (p1 == Glu && p2 == Gln)) ||
					((pep1[i] == Lys && p2 == Glu) || (p1 == Glu && p2 == Lys)) ||
					((pep1[i] == Ile && p2 == Asn) || (p1 == Asn && p2 == Ile)) ||
					((pep1[i] == Leu && p2 == Asn) || (p1 == Asn && p2 == Leu)))
				{
					dxy = 0.5;
				}
				else if (ox_met_aa > 0 && ((p1 == ox_met_aa && p2 == Phe) ||
					(p1 == Phe  && p2 == ox_met_aa)))
				{
					dxy = 0;
				}
			}

			v1 = old_row[j] + dxy;
			v2 = new_row[j - 1] + 1;
			v3 = (p2_pos<pep_len2 && j<2 * max_width) ? old_row[j + 1] + 1 : 9999;

			new_row[j] = v1;
			if (new_row[j]>v2)
				new_row[j] = v2;
			if (new_row[j]>v3)
				new_row[j] = v3;
		}

		float *tmp;
		tmp = old_row;
		old_row = new_row;
		new_row = tmp;
	}

	return old_row[pep_len2 - pep_len1 + max_width];
}





// changes all the amino acids to their original form (without PTMs)
void Peptide::convert_to_org(const Config *config)
{
	const vector<int>& org_aa = config->get_org_aa();
	unsigned int i;
	for (i = 0; i < amino_acids.size(); i++)
		amino_acids[i] = org_aa[amino_acids[i]];
}










