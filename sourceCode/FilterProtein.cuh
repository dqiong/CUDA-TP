
#include "includes.h"
#include "MonoSpectrum.h"

using namespace std;

double findKey(double b, map<double, int>& m);
void calWeight(const vector<double>& peaks, const vector<double>& masses);
void exactConvolution(const vector<double>& a, const vector<double>& b, double e);
void approxConvolution(const vector<double>& a, const vector<double>& b, double e);
int approxConvolution2(const vector<double>& masses, const vector<double>& peaks, double e);
vector<int> restrictedConvolution(vector<double> a, vector<double> b, double p, double delta, double e);
int msFilter(const vector<double>& a, const vector<double>& b, double e);
void msFilterOneSpectra(const vector<vector<double>>& a, const vector<double>& b, double e);

//void filterProtein(vector<vector<mass_t>>& cutMasses, vector<MonoSpectrum>& spectrum, vector<vector<mass_t>>& cutPeaksMonoMass);
//__device__ bool findVal(int* a, int value,int m);
//__global__ void msFiltergpu(int* cutMass_begin, int* cutMass_end, double* a,double* b,int m);
//__global__ void testt();

void testTime(vector<vector<mass_t>>& cutMasses, vector<vector<mass_t>>& cutPeaksMonoMass);