#include "stdafx.h"
#include <iostream>
#include "Cutpoint.h"


/**
* Constructor. Sets up internal tables. Calculates cutpoints.
*
* @param f tabulated d.f
* @param x arrays of x
* @param m number of cutpoints
*/
void Cutpoint::Initialize(std::vector<double>& x, std::vector<double>& f, int m)
{
	_X = x;
	_F = f;

	_m = m;
	_I.resize(m);

	int j = -1;
	double a = 0;

	for (int i = 0; i<_m; i++)
	{
		while (a <= i && ++j < _F.size())
		{
			a = _m*_F[j];
		}
		_I[i] = j;
	}
}

int Cutpoint::return_I(int index)
{
	return _I[index];
}

int Cutpoint::return_z()
{
	return z;
}