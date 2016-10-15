
#include "Config.h"
#include "auxfun.h"



void Config::init_with_defaults()
{
	unsigned int i;

	mass_spec_type = ESI_MASS_SPEC; // default type
	resource_dir = "";
	config_file = "";

	tolerance = 0.5;
	pm_tolerance = 0.5;

	local_window_size = 200;
	max_number_peaks_per_local_window = 40;
	number_of_strong_peaks_per_local_window = 12;

	resource_dir = "Models";

	init_standard_aas();

	session_aas = standard_aas;

	// these all operate on the original aas
	session_tables.init_for_standard_aas();

	// insert labels of original aas
	label2aa.clear();
	const vector<string>& aa2label = get_aa2label();
	for (i = 0; i < aa2label.size(); i++)
		label2aa.insert(STRING2INT_MAP::value_type(aa2label[i], i));



	string PTM_file = "TopDownPTMs.txt";

	ind_read_PTM_file = false;

	read_PTM_file(PTM_file);
}

void Config::init_standard_aas()
{
	int i;
	standard_aas.clear();
	for (i = Ala; i <= Val; i++)
		standard_aas.push_back(i);
}






/***********************************************************
// returns the idx of an aa from its label
// -1 if label is not found
***********************************************************/
int Config::get_aa_from_label(const string& label) const
{
	STRING2INT_MAP::const_iterator iter = label2aa.find(label);

	if (iter == label2aa.end())
		return -1;

	return (*iter).second;
}




bool read_mass_type_line(const char* prefix, char *line, mass_t& val)
{
	int len = strlen(prefix);
	if (strncmp(prefix, line, len))
		return false;
	istringstream is(line + len);
	is >> val;
	return true;
}




/************************************************************
outputs the selected aas for the different regions
*************************************************************/
void Config::print_table_aas(const ConversionTables& table,
	const vector<int>& aas) const
{
	unsigned int i;
	cout << aas.size() << " amino acids:" << endl;

	for (i = 0; i < aas.size(); i++)
		cout << setw(6) << left << aas[i] << setprecision(4) << right << fixed << setw(8)
		<< table.get_aa2mass(aas[i]) << "   " << left <<
		table.get_aa2label(aas[i]) << endl;
}

void Config::print_session_aas() const
{
	cout << endl << "AMINO ACIDS" << endl;
	print_table_aas(session_tables, session_aas);
	cout << "N_TERM " << session_tables.get_aa2mass(N_TERM) << endl;
	cout << "C_TERM " << session_tables.get_aa2mass(C_TERM) << endl;
}






