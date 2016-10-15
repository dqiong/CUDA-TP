

#ifndef __AUXFUN_H__
#define __AUXFUN_H__

#include "includes.h"


void end(int status = 0);

/****************************************************************
Gets a list of k sets (each set is a list of ints)
return k lists of length n (for position i this is the cross
product i) for the n posisble cross products
*****************************************************************/
void generateAllCrossProducts(const vector< int >& sizesList,
	vector< vector<int> >& products);

mass_t ppm_val(mass_t offset, mass_t total_mass);


void rand_seed(unsigned int init = 0);



/* Returns random uniform number */
double my_random();


void choose_k_from_n(int k, int n, vector<int>& idxs);

void generate_all_permutations(const vector<int>& org_vector,
	vector< vector<int> >& permutations);

// returns the minimal x for which the cumulative probability
// P(X<x)>= target_prob, assuming X~bin(n,p)
int get_min_number_from_binomial_prob(int n, double p, double target_prob);

//创建直方图函数
template<class T>
void create_histogram(vector<T>& vals, int num_bins, T min_val,
	T max_val, ostream& os = cout)
{
	T bin_size = (max_val - min_val) / (T)(num_bins);
	T one_over_bin = 1.0 / bin_size;
	int i;
	vector<int> counts;
	vector<float> percents;
	counts.resize(num_bins, 0);

	for (i = 0; i<vals.size(); i++)
	{
		if (vals[i]<min_val)
		{
			counts[0]++;
			continue;
		}

		int bin_idx = num_bins - 1;
		if (vals[i]<max_val)
			bin_idx = (int)(one_over_bin*(vals[i] - min_val));

		counts[bin_idx]++;
	}

	T v = min_val;
	int tc = 0;
	for (i = 0; i<num_bins; i++)
	{
		os << setw(4) << setprecision(2) << left << v << " - ";
		v += bin_size;
		os << setw(4) << left << v << "  " << setw(6) << right << counts[i] << "  ";
		os << setw(6) << left << setprecision(4) << (float)counts[i] / (float)vals.size() << endl;
		tc += counts[i];
	}

	os << "Total:       " << setw(6) << right << tc << "  " << setw(6) << setprecision(4) << left << (float)tc / (float)vals.size() << endl;
}


template<class T>
void create_histogram(const vector<T>& vals, const vector<T>& separator_vals,
	vector<int>& counts, ostream& os = cout)
{

	int i;

	vector<float> percents;
	counts.resize(separator_vals.size() + 1, 0);

	for (i = 0; i<vals.size(); i++)
	{
		int j;
		for (j = 0; j<separator_vals.size(); j++)
		if (vals[i] <= separator_vals[j])
			break;

		counts[j]++;
	}

	/*	cout << "VALS: " << vals.size() << endl;


	int tc=0;

	for (i=0; i<counts.size(); i++)
	{
	T v = (i==0) ? 0 : separator_vals[i-1];
	os << setw(4) <<  setprecision(2) << left << v << " - ";
	v = ( i == counts.size() -1 ) ? 0 : separator_vals[i];
	os <<  setw(4) << left << v  << "  " <<  setw(6) << right << counts[i] << "  ";
	os << setw(6) << left << setprecision(4) << (float)counts[i]/(float)vals.size() << endl;
	tc+= counts[i];
	}

	os << "Total:       " << setw(6) << right << tc << "  " << setw(6) << setprecision(4) << left << (float)tc/(float)vals.size() << endl;
	*/
}


//计算平均值和方差
template<class T>
void calc_mean_sd(const vector<T>& v, T *mean, T *sd)
{
	T m = 0, var = 0;
	unsigned int i;

	if (v.size() == 0)
	{
		*mean = 0;
		*sd = 0;
		return;
	}

	if (v.size() == 1)
	{
		*mean = v[0];
		*sd = 0;
	}

	for (i = 0; i<v.size(); i++)
		m += v[i];

	m /= v.size();

	for (i = 0; i<v.size(); i++)
		var += (v[i] - m)*(v[i] - m);

	var /= v.size();

	*mean = m;
	*sd = sqrt(var);
}

//把a和b中的相同元素去除,a中剩下的赋给diff
template<class T>
void set_minus(const vector<T>& a, const vector<T>& b, vector<T>& diff)
{
	int i;
	diff.clear();
	for (i = 0; i<a.size(); i++)
	{
		int j;
		for (j = 0; j<b.size(); j++)
		if (a[i] == b[j])
			break;
		if (j == b.size())
			diff.push_back(a[i]);
	}
}

//把a和b中的相同的元素赋给overlap
template<class T>
void set_overlap(const vector<T>& a, const vector<T>& b, vector<T>& overlap)
{
	int i;
	overlap.clear();
	for (i = 0; i<a.size(); i++)
	{
		int j;
		for (j = 0; j<b.size(); j++)
		if (a[i] == b[j])
			break;
		if (j<b.size())
			overlap.push_back(a[i]);
	}
}


struct PermutePair {
	bool operator< (const PermutePair& other) const
	{
		return randVal<other.randVal;
	}
	int    orgIdx;
	double randVal;
};
//把vec中的数随机打乱
template<class T>
void permuteVector(vector<T>& vec)
{
	vector<PermutePair> idxPairs;
	idxPairs.resize(vec.size());
	unsigned int i;
	for (i = 0; i<vec.size(); i++)
	{
		idxPairs[i].orgIdx = i;
		idxPairs[i].randVal = my_random();
	}

	sort(idxPairs.begin(), idxPairs.end());
	vector<T> new_vec;
	new_vec.resize(vec.size());
	for (i = 0; i<vec.size(); i++)
		new_vec[i] = vec[idxPairs[i].orgIdx];

	vec = new_vec;
}

#endif 
