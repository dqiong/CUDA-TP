
#ifndef __MONOFRAGMENTS_H__
#define __MONOFRAGMENTS_H__

#include "includes.h"

typedef enum FragmentDir  { PREFIX, SUFFIX } FragmentDir;
typedef enum FragmentType { STRONG, WEAK } FragmentType;

class MonoFragment {
	friend class FragmentSet;
public:
	MonoFragment() : label(""), direction(PREFIX), type(STRONG), offset(0) {};
	MonoFragment(string _label, FragmentDir _direction, FragmentType _type, mass_t _offset) :
		label(_label), direction(_direction), type(_type), offset(_offset) {};

	// totalMonoMassWith18
	mass_t calcExpFragMass(mass_t cutMass, mass_t totalMonoMassWith18) const
	{
		mass_t actualCutMass = (direction == PREFIX) ? cutMass : totalMonoMassWith18 - cutMass;
		return actualCutMass + offset;
	}

	mass_t calcExpCutMass(mass_t fragMass, mass_t totalMonoMassWith18) const
	{
		if (direction == PREFIX)
			return fragMass - offset;

		return (totalMonoMassWith18 - fragMass) + offset;
	}

	FragmentDir getDirection() const { return direction; }

	void print(ostream& os = cout);

	string getLabel() const { return label; }

private:
	string		 label;
	FragmentDir  direction;
	FragmentType type;
	mass_t		 offset;

};


class FragmentSet {
public:
	void calcFragmentSetMonoMasses(mass_t cuttMass, mass_t totalMonoMass, vector<mass_t>& masses) const;

	void initECD();

	int	getNumFragments() const { return fragments.size(); }

	const MonoFragment& getFragment(int idx) const { return fragments[idx]; }

	void printFragments() const;

	void printFragmentNames(ostream& os = cout) const;

	const vector<int>& getStrongFragmentIdxs() const { return strongFragmentIdxs; }

private:
	vector<MonoFragment> fragments;

	vector<int>	strongFragmentIdxs; // idxs of fragments to be used to extrapolate prefix masses

};

#endif