
#include "MonoFragments.h"




void FragmentSet::initECD()
{
	strongFragmentIdxs.clear();
	strongFragmentIdxs.push_back(0);
	strongFragmentIdxs.push_back(1);
	fragments.push_back(MonoFragment("c", PREFIX, STRONG, 17.02635));
	fragments.push_back(MonoFragment("z", SUFFIX, STRONG, -16.01872)); // includes 18 Da terminal

	fragments.push_back(MonoFragment("c-H2O", PREFIX, WEAK, -1.07928));
	fragments.push_back(MonoFragment("c-NH3", PREFIX, WEAK, -0.000197));
	fragments.push_back(MonoFragment("z-H2O", SUFFIX, WEAK, -34.029283)); // includes 18 Da terminal
	fragments.push_back(MonoFragment("z-NH3", SUFFIX, WEAK, -33.045267)); // includes 18 Da terminal

	//	fragments.push_back(MonoFragment("b",PREFIX,WEAK,    0));
	//	fragments.push_back(MonoFragment("y",SUFFIX,WEAK,    0)); // includes 18 Da terminal
	//	fragments.push_back(MonoFragment("c-2H2O",PREFIX,WEAK, -19.089843));
	//	fragments.push_back(MonoFragment("z-2H2O",SUFFIX,WEAK,  -52.039846)); // includes 18 Da terminal
	//	
}



void FragmentSet::calcFragmentSetMonoMasses(mass_t cutMass, mass_t totalMonoMass,
	vector<mass_t>& fragMasses) const
{
	fragMasses.clear();
	vector<MonoFragment>::const_iterator it;
	for (it = fragments.begin(); it != fragments.end(); it++)
		fragMasses.push_back(it->calcExpFragMass(cutMass, totalMonoMass));

}



void FragmentSet::printFragmentNames(ostream& os) const
{
	unsigned int i;
	for (i = 0; i < fragments.size() - 1; i++)
		os << fragments[i].label << "\t";
	os << fragments[i].label << endl;
}
