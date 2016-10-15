
#include "Config.h"
#include "auxfun.h"

//实现上面两个文件中的函数
/**************************************************************
Reads a list of PTMS from a file
C    57.021464   FIXED   ALL    C+57    Carbamidomethyl
***************************************************************/
void Config::read_PTM_file(const string& file)
{
	if (ind_read_PTM_file)
		return;

	char buff[256];
	FILE *stream = fopen(file.c_str(), "r");   //标记一下
	if (!stream)
	{
		cout << "Error: couldn't open PTM file: " << file << endl;
		end(1);
	}

	all_PTMs.clear();

	while (fgets(buff, 256, stream))
	{
		char aa_labels[64];
		char offset_str[64];
		char type_str[64];
		char region_str[64];
		char symbol[64];
		char name[64];

		int num_args_read = sscanf(buff, "%s %s %s %s %s %s", aa_labels, offset_str, type_str, region_str, symbol, name);
		if (num_args_read != 6)
			continue;

		if (aa_labels[0] == '#')
			continue;

		int i;
		int len = strlen(aa_labels);
		for (i = 0; i < len; i++)
		{
			if (aa_labels[i] == ':' || aa_labels[i] == ',')
				aa_labels[i] = ' ';
		}

		vector<int> ptmAminoAcids;
		istringstream iss(aa_labels);
		while (1)
		{
			int aa = -1;
			string single_label = "";
			iss >> single_label;

			if (single_label.length() < 1)
				break;

			aa = get_aa_from_label(single_label.c_str());

			if (!strcmp("N_TERM", single_label.c_str()))
				aa = N_TERM;

			if (!strcmp("C_TERM", single_label.c_str()))
				aa = C_TERM;

			if (aa < 0)
			{
				cout << "Unkown aa " << single_label << "  in PTM:" << endl << buff << endl;
				end(1);
			}
			else
				ptmAminoAcids.push_back(aa);
		}

		mass_t offset = (mass_t)atof(offset_str);
		int type = (!strcmp(type_str, "OPTIONAL")) ? PTM_OPTIONAL : PTM_FIXED;

		int region = -1;
		int position = 0;
		if (!strcmp("ALL", region_str))
		{
			region = PTM_ALL;
		}
		else if (!strcmp("C_TERM", region_str))
		{
			region = PTM_C_TERMINAL;
		}
		else if (!strcmp("N_TERM", region_str))
		{
			region = PTM_N_TERMINAL;
		}
		else
		{
			position = atoi(region_str);
			if (position<10 && position>-10)
				region = PTM_POSITION;
		}

		if (region < 0)
		{
			cout << "Error: bad PTM region : " << region_str << endl;
			end(1);
		}

		// add PTMs
		PTM newPtm;

		newPtm.applicableAAs = ptmAminoAcids;
		newPtm.delta = offset;
		newPtm.type = type;
		newPtm.region = region;
		newPtm.position = position;
		newPtm.label = symbol;
		newPtm.name = name;

		all_PTMs.add_PTM(newPtm);

	}

	// sanity check (no two PTMs can have the same label)
	vector<string> used_labels;

	int i;
	for (i = 0; i < all_PTMs.get_num_PTMs(); i++)
	{
		unsigned int j;
		for (j = 0; j < used_labels.size(); j++)
		{
			if (used_labels[j] == all_PTMs.list[i].label)
			{
				cout << "Error: " << used_labels[j] << " is already used!" << endl;
				end(1);
			}
		}
		used_labels.push_back(all_PTMs.list[i].label);
	}


	// update the conversion vectors according to the PTMS

	session_tables.init_for_standard_aas();


	// create a mapping between label and aa
	const vector<string>& aa2label = session_tables.get_aa2label();
	const int size = aa2label.size();
	for (i = 0; i < size; i++)
	{
		label2aa.insert(STRING2INT_MAP::value_type(aa2label[i], i));
	}


	ind_read_PTM_file = true;

	//	this->print_supported_PTMs();
}



/**************************************************************
// returns the idx of the label
// -1 is returned if label is not found
***************************************************************/
int PTM_list::get_PTM_idx(const string& label) const
{
	unsigned int i;

	for (i = 0; i < list.size(); i++)
	{
		const string & l_label = list[i].label;
		if (list[i].label == label)
			return i;
	}

	return -1;
}


/**************************************************************
// adds a PTM, makes sure that it wasn't already added
// (if already added, returns false
***************************************************************/
bool PTM_list::add_PTM(const PTM& ptm)
{
	unsigned int i;
	for (i = 0; i < list.size(); i++)
	{
		if (ptm.label == list[i].label)
			return false;
	}

	if (ptm.label.length() < 2)
	{
		//		cout << "Error: PTM label must be longer (form should be X-dd or X+dd)\n";
		//		end(1);
	}
	list.push_back(ptm);
	return true;
}







/********************************************************************
Gets a text line with all the PTM labels that are to be used in this
run. Updates the conversion tables and PTM lists appropriately
Only fixed mods cause a change to the masses (the optional mod are
assumed to already have the correction).
This function erases the effects of any previous PTMs selected!
*********************************************************************/
void Config::apply_fixed_PTMs(char *ptm_line)
{
	// replace all colons with white space
	int i;
	int ptm_line_length = strlen(ptm_line);
	for (i = 0; i < ptm_line_length; i++)
	if (ptm_line[i] == ':')
		ptm_line[i] = ' ';

	istringstream ptm_stream(ptm_line);
	string ptm_label;

	// create lists of AAs for different regions
	session_aas = standard_aas;

	int num_fixed_C_TERM_PTMs = 0, num_fixed_N_TERM_PTMs = 0;

	while (ptm_stream >> ptm_label)
	{
		int idx;

		idx = all_PTMs.get_PTM_idx(ptm_label);
		if (idx < 0 || all_PTMs.list[idx].type != PTM_FIXED)
		{
			continue;
		}

		if (idx >= 0)
		{
			const PTM& ptm = all_PTMs.list[idx];

			unsigned int j;
			for (j = 0; j < selected_fixed_PTM_idxs.size(); j++)
			{
				if (selected_fixed_PTM_idxs[j] == idx)
					break;
			}

			if (j < selected_fixed_PTM_idxs.size())
				continue;

			selected_fixed_PTM_idxs.push_back(idx);

			if (ptm.region == PTM_ALL)
			{
				session_tables.make_fixed_mod(ptm.applicableAAs[0], ptm.delta);
				continue;
			}

			cout << "ERROR: Fixed PTMs must be in region PTM_ALL!" << endl;
			end(1);
		}



		cout << "Error: No support for PTM: " << ptm_label << endl;
		end(1);
	}

	sort(session_aas.begin(), session_aas.end());

	//	calc_aa_combo_masses(2);

	//	set_aa_variants();
}


ostream& operator << (ostream& os, const PTM& ptm)
{
	os << setw(8) << ptm.label << " " << setw(12) << fixed << right << ptm.delta << "  ";
	os << setw(12) << left;
	if (ptm.region == PTM_ALL)
	{
		os << "ALL";
	}
	else if (ptm.region == PTM_N_TERMINAL)
	{
		os << "N_TERM";
	}
	else if (ptm.region == PTM_C_TERMINAL)
	{
		os << "C_TERM";
	}
	else if (ptm.region == PTM_POSITION)
	{
		os << ptm.position;
	}

	os << setw(10) << left << ((ptm.type == PTM_OPTIONAL) ? "OPTIONAL " : "FIXED ");

	if (!ptm.name.empty())
		os << " (" << ptm.name << ")";

	ConversionTables conv;
	conv.init_for_standard_aas();
	os << "\t";
	for (int i = 0; i < ptm.applicableAAs.size(); i++){  //输出对应的氨基酸
		int k = ptm.applicableAAs[i];
		os << conv.get_aa2char()[k] << ",";
	}
	//delete *config;

	return os;
}


void Config::print_supported_PTMs() const
{
	unsigned int i;

	cout << endl << "MODIFICATIONS ( " << all_PTMs.list.size() << " )" << endl;
	for (i = 0; i < all_PTMs.list.size(); i++)
		cout << left << setw(3) << i + 1 << " " << all_PTMs.list[i] << endl;

	cout << endl;
}

