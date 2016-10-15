/***************************************************************************
* Title:          Config.h
* Author:         Ari Frank
* Copyright (c) 2009 The Regents of the University of California
* All Rights Reserved
* See file LICENSE for details.
***************************************************************************/

#ifndef __CONFIG_H__
#define __CONFIG_H__

#include "ConversionTables.h"
#include "BasicDataStructs.h"
#include "includes.h"

typedef enum PTM_REGIONS { PTM_ALL, PTM_N_TERMINAL, PTM_C_TERMINAL, PTM_POSITION } PTM_REGIONS;

typedef enum PTM_TYPES { PTM_FIXED, PTM_OPTIONAL } PTM_TYPES;


// the Config class holds all configuration variables that are used
// by the models and programs: aa's, PTMs, thresholds etc


// label for fixed PTM stays the same as the orginal amino acid
struct PTM {
	PTM() : delta(0) {};

	vector<int> applicableAAs; // holds the amino acids that can hold this PTM

	mass_t delta; // relative to unmodified amino acid

	int type;   // PTM_FIXED, PTM_OPTIONAL

	int region; // PTM_ALL , PTM_N_TERMINAL, PTM_C_TERMINAL, PTM_POSITION

	int position; // used only for specific position like Q-17 at the +1 position

	string label;

	string name;

	bool wasApplied; // for fixed
};

ostream& operator << (ostream& os, const PTM& ptm);



struct PTM_list {
	friend class Config;
public:

	void clear() { list.clear(); }

	int get_num_PTMs() const { return list.size(); }

	const vector<PTM>& getAllPTMs() const { return list; }


	// returns the idx of the label
	// -1 is returned if label is not found
	int get_PTM_idx(const string& label) const;

private:

	// adds a PTM, makes sure that it wasn't already added
	// (if already added, returns false
	bool add_PTM(const PTM& ptm);

	vector<PTM> list;
};


// A class that holds a map of mass ranges
class MassRangeMap {
public:
	MassRangeMap() : was_initialized(0) { clear(); }

	void clear(mass_t max_mass = 400);

	// adds a mass range to the map, if already covered, does nothinh
	void add_range(mass_t min_mass, mass_t range);

	// adds a new set of ranges *shift_size* away from the current existing ranges
	void add_shifted_ranges(mass_t shift_size);

	// checks if the given mass is in one of the maps ranges
	bool is_covered(mass_t mass) const
	{
		if (mass > max_map_mass)
			return true;

		MASS_T_MAP::const_iterator it;
		it = ranges.lower_bound(mass);
		if (it != ranges.end() && it != ranges.begin())
		{
			if (it->first == mass)
				return true;

			it--;
			if (it->first <= mass && it->second >= mass)
				return true;
		}
		return false;
	}

	void print(ostream& os = cout) const;

	void read_ranges(istream& is);
	void write_ranges(ostream& os = cout) const;

	void set_was_initialized(int flag) { was_initialized = flag; }
	int  get_was_initialized() const { return was_initialized; }

private:
	int was_initialized;
	mass_t max_map_mass;
	MASS_T_MAP ranges;
};




class Config {
public:

	// sets the values for all defined PTMs
	//从文件中读取ptm相关信息
	// any additional PTMs can only be user defined


	void init_with_defaults();

	void init_regional_fragment_set_defaults(int type_set = 0, int max_charge = 5);

	void print_supported_PTMs() const;

	void apply_fixed_PTMs(char *ptm_line);

	void read_PTM_file(const string& file);

	int get_num_PTMs() const { return all_PTMs.list.size(); }


	// returns the idx of the label
	// -1 is returned if label is not found
	int get_PTM_idx(const string& label) const { return all_PTMs.get_PTM_idx(label); }
	const PTM& get_PTM(int idx) const { return all_PTMs.list[idx]; }
	const vector<PTM>& getAllPTMs() const { return all_PTMs.getAllPTMs(); }

	void fill_allowed_double_edges(bool allow_all = false);

	// returns the idx of an aa from its label
	// -1 if label is not found
	int get_aa_from_label(const string& label) const;


	const ConversionTables& get_session_tables() const { return session_tables; }

	const vector<int>& get_session_aas() const { return session_aas; }
	const vector<int>& get_char2aa() const { return session_tables.get_char2aa(); }
	const vector<string>& get_aa2label() const { return session_tables.get_aa2label(); }
	const vector<mass_t>& get_aa2mass() const { return session_tables.get_aa2mass(); }
	const vector<mass_t>& get_char2mass() const { return session_tables.get_char2mass(); }
	const vector<int>&    get_org_aa() const { return session_tables.get_org_aa(); }

	int    get_digest_type() const { return digest_type; }
	void   set_digest_type(int type);

	const vector<int>& get_n_term_digest_aas() const { return n_term_digest_aas; }
	const vector<int>& get_c_term_digest_aas() const { return c_term_digest_aas; }
	int   get_num_n_term_digest_aas() const { return n_term_digest_aas.size(); }
	int   get_num_c_term_digest_aas() const { return c_term_digest_aas.size(); }
	bool  is_n_digest_aa(int aa) const { unsigned int i; for (i = 0; i < n_term_digest_aas.size(); i++) if (aa == n_term_digest_aas[i]) break; return (i < n_term_digest_aas.size()); }
	bool  is_c_digest_aa(int aa) const { unsigned int i; for (i = 0; i < c_term_digest_aas.size(); i++) if (aa == c_term_digest_aas[i]) break; return (i < c_term_digest_aas.size()); }


	int    get_need_to_estimate_pm() const { return need_to_estimate_pm; }
	void   set_need_to_estimate_pm(int val) { need_to_estimate_pm = val; }
	int    get_need_to_normalize() const { return need_to_normalize; }
	void   set_need_to_normalize(int val) { need_to_normalize = val; }
	int    get_mass_spec_type() const { return mass_spec_type; }

	mass_t get_tolerance() const { return tolerance; }
	mass_t get_pm_tolerance() const { return pm_tolerance; }

	void set_tolerance(mass_t t) { tolerance = t; }
	void set_pm_tolerance(mass_t t) { pm_tolerance = t; }
	void set_tolerances(mass_t t) { tolerance = t; pm_tolerance = t; }

	mass_t get_local_window_size() const { return local_window_size; }
	int    get_max_number_peaks_per_local_window() const { return max_number_peaks_per_local_window; }
	void   set_max_number_peaks_per_local_window(int n) { max_number_peaks_per_local_window = n; }
	int    get_number_of_strong_peaks_per_local_window() const { return number_of_strong_peaks_per_local_window; }
	string get_resource_dir() const { return resource_dir; }
	string get_config_file() const { return config_file; }
	void set_resource_dir(string& _resource_dir) { resource_dir = _resource_dir; }
	void set_config_file(string& _config_file) { config_file = _config_file; }


	void print_session_aas() const;
	void print_fragments(ostream &os) const;
	void print_regional_fragment_sets(ostream& os = cout) const;
	void read_fragments(istream& is);
	void read_regional_fragment_sets(istream& is);
	void clone_regional_fragment_sets(int source_charge, int target_charge);


private:
	//	ConversionTables original_tables;  // with default values

	// These conversion tables represent the parameters after PTM modifications
	// All tables have the same aa indices however the actual masses might be
	// different (due to terminal additions for instance)
	ConversionTables session_tables;

	vector<int> standard_aas; // the 20 unmodified amino acids
	vector<int> session_aas; //  all aas that can be part of a peptide (including terminal mods)
	// all the variants that can be used for each amino acid (e.g. M-> M,M+16,M+32)
	// maps all labels to their aa_idx
	STRING2INT_MAP label2aa;



	// MASS SPEC TYPE
	int		   mass_spec_type; // type of machine, might influence the way things are done
	// this parameter should be part of the model file, and not changed.

	int		    digest_type; // 0 - nothing , 1 Trypsin

	vector<int> n_term_digest_aas;
	vector<int> c_term_digest_aas;

	int        need_to_estimate_pm; // 0 - used file pm, 1 - use original pm

	int		   need_to_normalize; // 0 - don't normalize intensities, 1 - do normalize (sum of intensities = m_over_z)

	string     resource_dir;  // path to direcotry where resource files can be found
	string     config_file;

	// PTM vectors (hold info on all supported PTMs)
	bool       ind_read_PTM_file;

	PTM_list    all_PTMs;    // region PTM_ALL
	vector<int> selected_fixed_PTM_idxs;

	// Tolerance variables
	mass_t tolerance;       // tolerance for consecutive peaks
	mass_t pm_tolerance;

	mass_t local_window_size;
	int	 max_number_peaks_per_local_window;
	int  number_of_strong_peaks_per_local_window;

	void init_standard_aas();

	void print_table_aas(const ConversionTables& table,
		const vector<int>& aas) const;



};
#endif


