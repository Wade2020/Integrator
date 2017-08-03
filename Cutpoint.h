#pragma once
#ifndef CUTPOINT_H
#define CUTPOINT_H

#include <cmath>
#include <vector>


/**
* \ingroup monte_carlo_grp
*
* \class Cutpoint
*
* \brief Generates random variable with tabulated distribution
* fucntion using Cutpoint method (Fishman G.S. Algorithm CM/CMSET)
*
*/
class Cutpoint
{
public:

	//! Constructor
	void Initialize(std::vector<double>& f, std::vector<double>& x, int m);

	//! return value of random variable with the tabulated d.f.
	double operator()(double ksi);
	double operator()(double ksi, std::vector<int> _I, std::vector<double> _X, std::vector<double> _F, int _m);

	int return_I(int index);
	int return_z();
	int co;

private:

	int                 _m;  //!< number of cutpoints
	std::vector<int>    _I;  //!< cutpoint index array

	std::vector<double> _F;  //!< tabulated d.f. \f$ P{\xi < x } = F(x) \f$
	std::vector<double> _X;	 //!< array of x

	int z;
};


/**
* returns value of random variable with d.f. \f$ F(x) \f$
*
* @param  ksi uniformly distributed on [0,1) random variate
*
* @return random variable with the tabulated d.f. \f$ F(x) \f$
*/

inline double Cutpoint::operator()(double ksi)
{
	z = _I[static_cast<int>(floor(_m*ksi))];
	co = 0;
	while (ksi > _F[z]) {
		z++;
		co++;
	}

	return _X[z];
}

inline double Cutpoint::operator()(double ksi, std::vector<int> _I, std::vector<double> _X, std::vector<double> _F, int _m)
{
	z = _I[static_cast<int>(floor(_m*ksi))];
	co = 0;
	while (ksi > _F[z]) {
		z++;
		co++;
	}

	return _X[z];
}



#endif