#include "auxfun.h"

void end(int status)
{
	system("pause");
	exit(status);
}

static unsigned int SEED;

void rand_seed(unsigned int init)   {
	if (init != 0)
	{
		SEED = init;
	}
	else
	{
		time_t ltime;
		unsigned int t = (unsigned int)time(&ltime);

		SEED = t;
	}
}


/* Returns random uniform number 返回一个随机数,是一个0.几的6位小数*/
double my_random()
{
	static unsigned int a = 1588635695, m = 4294967291U, q = 2, r = 1117695901;

	SEED = a*(SEED % q) - r*(SEED / q);
	return ((double)SEED / (double)m);
}



// chooses k numbers from 0,...,n-1 (unique)
void choose_k_from_n(int k, int n, vector<int>& idxs)
{
	int i;
	idxs.clear();
	idxs.reserve(k);

	if (k>n)
	{
		cout << "Error: choose " << k << " from " << n << " !" << endl;
		end(1);
	}

	i = 0;
	while (i<k)
	{
		int idx = (int)(my_random() * n);
		int j;
		for (j = 0; j<i; j++)
		if (idx == idxs[j])
			break;

		if (j<i)
			continue;

		idxs.push_back(idx);
		i++;
	}

	sort(idxs.begin(), idxs.end());
}

mass_t ppm_val(mass_t offset, mass_t total_mass)
{
	return (offset / total_mass) * 1000000;
}



/*************************************************************
finds all the permutaitons of n elements, repeated elements
are allowed and do not create redundant permutations.
产生n个元素的所有全排列,允许重复元素,但是不产生冗余的排列
**************************************************************/
void generate_all_permutations(const vector<int>& org_vector,
	vector< vector<int> >& permutations)
{
	unsigned int i;
	vector<int> counts, symbols;
	permutations.clear();

	if (org_vector.size() == 0)
		return;

	counts.clear();
	symbols.clear();

	// create vector with symbols and their counts
	symbols.push_back(org_vector[0]);
	counts.push_back(1);

	for (i = 1; i<org_vector.size(); i++)
	{
		unsigned int j;
		for (j = 0; j<counts.size(); j++)
		{
			if (org_vector[i] == symbols[j])
			{
				counts[j]++;
				break;
			}
		}

		if (j == counts.size())
		{
			symbols.push_back(org_vector[i]);
			counts.push_back(1);
		}
	}

	vector<int> next_sym_idx, perm;
	int n = org_vector.size(); // total number of elements
	int k = counts.size(); // total number of element types
	next_sym_idx.resize(n, 0);
	perm.resize(n, -1);
	int d = 0;

	while (1)
	{
		while (next_sym_idx[d]<k && counts[next_sym_idx[d]] == 0)
			next_sym_idx[d]++;

		if (next_sym_idx[0] == k)
			break;

		if (next_sym_idx[d] >= k)
		{
			next_sym_idx[d] = 0;
			d--;
			counts[next_sym_idx[d]]++;
			next_sym_idx[d]++;
			continue;
		}

		// add symbol
		perm[d] = symbols[next_sym_idx[d]];
		counts[next_sym_idx[d]]--;
		d++;

		if (d == n)
		{
			permutations.push_back(perm);
			//		int k;
			//		for (k=0; k<perm.size(); k++)
			//			cout << perm[k] << " ";
			//		cout << endl;

			d--;
			counts[next_sym_idx[d]]++;
			next_sym_idx[d]++;
		}
	}
}

// returns the minimal x for which the cumulative probability
// P(X<x)>= target_prob, assuming X~bin(n,p)
//返回最小的X的累积概率,X服从二项分布
int get_min_number_from_binomial_prob(int n, double p, double target_prob)
{
	const double one_minus_p = 1.0 - p;

	double pow_p = 1.0;
	double pow_1_minus_p = pow(one_minus_p, n);

	double sum_prob = pow_1_minus_p;
	double bin_coef = 1.0;
	double pow_val = pow_1_minus_p;

	//	cout << 0 << " " <<  pow_val << " " << pow_val << endl;
	int b = 0;
	while (sum_prob<target_prob)
	{
		b++;
		bin_coef *= (double)(n - b + 1);
		bin_coef /= (double)b;

		pow_val *= p;
		pow_val /= one_minus_p;

		double prob = bin_coef * pow_val;
		sum_prob += prob;

		//	cout << b << " " << prob << " " << sum_prob << endl;
	}
	return b;
}

/****************************************************************
Gets a list of k sets (each set is a list of ints)
return k lists of length n (for position i this is the cross
product i) for the n posisble cross products
生成所有排列,如第一位是2,那么这个位置上只能0或者1,貌似没啥用
*****************************************************************/
void generateAllCrossProducts(const vector<int>& sizesList,
	vector< vector<int> >& crossProducts)
{
	const int numSizes = sizesList.size();
	const int lastSize = sizesList.size() - 1;
	vector<int> idxs;
	vector<int> v;

	v.resize(numSizes, 0);
	idxs.resize(numSizes, 0);
	crossProducts.clear();

	int d = 0;
	while (1)
	{
		if (d<0)
			break;

		if (d == lastSize)
		{
			int i;
			for (i = 0; i<sizesList[lastSize]; i++)
			{
				v[lastSize] = i;
				crossProducts.push_back(v);
			}
			d--;
			continue;
		}

		if (idxs[d] == sizesList[d])
		{
			idxs[d] = 0;
			d--;
			continue;
		}

		v[d] = idxs[d]++;
		d++;
	}

	/*	size_t i;
	for (i=0; i<sizesList.size(); i++)
	cout << sizesList[i] << " ";
	cout << endl;

	for (i=0; i<crossProducts.size(); i++)
	{
	size_t j;
	cout << i << " : ";
	for (j=0; j<crossProducts[i].size(); j++)
	cout << " " << crossProducts[i][j];
	cout << endl;
	}*/
}
