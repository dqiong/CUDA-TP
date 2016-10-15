

#ifndef __INCLUDES_H__
#define __INCLUDES_H__

#define _CRT_SECURE_NO_DEPRECATE 

#pragma warning (disable:4786)
#pragma warning (disable:4305)

#include <iostream>
#include <sstream>
#include <string>
#include <iomanip>
#include <fstream>
#include <algorithm>
#include <vector>
#include <cmath>
#include <fstream>
#include <map>
#include <time.h>
#include <iomanip>




using namespace std;
using std::vector;



// Improved vector that checks access violations - has overhead use only in DEBUG
#if defined(WIN32) && defined(_DEBUG)
template<class T, class A = allocator<T> >
class my_vector : public vector<T, A>
{
public:
	explicit my_vector(const A& al = A()) : vector<T, A>(al) {}
	explicit my_vector(size_type n, const T& v = T(), const A& al = A()) : vector<T, A>(n, v, al) {}

	my_vector(const my_vector& x) : vector<T, A>(x) {}
	my_vector(const_iterator first, const_iterator last, const A& al = A()) : vector<T, A>(first, last, al) {}

	const_reference operator[](size_type pos) const
	{
		return at(pos);
	}

	reference operator[](size_type pos)
	{
		return at(pos);
	}
};
#define vector my_vector
#endif //defined(_DEBUG) && defined(_WIN32)



#define NEG_INF -1E3

#define NUM_SIG_DIGITS 3

#define ESI_MASS_SPEC 1

typedef enum AminoAcids {
	N_TERM, C_TERM, Gap, Xle, Ala, Arg, Asn, Asp, Cys, Gln, Glu, Gly, His,
	Ile, Leu, Lys, Met, Phe, Pro, Ser, Thr, Trp, Tyr, Val
} AminoAcids;

#define MASS_H2O 18.010563
#define MASS_NH3 17.026547
#define MASS_NEUTRON 1.0086649101 


// Data types for common variables

typedef double mass_t;
typedef float  intensity_t;
typedef float  score_t;

typedef map< string, int, less<string> > STRING2INT_MAP;
typedef map< mass_t, mass_t, less<mass_t> > MASS_T_MAP;
typedef map< int, int, less<int> > INT_MAP;
typedef map< mass_t, int, less<mass_t> > MASS_T2INT_MAP;


#endif


