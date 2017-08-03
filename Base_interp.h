#pragma once
#include <vector>
#include <algorithm>

struct Base_interp
{
	// Abstract base class used by all interpolation routines in this chapter. Only the routine interp is called directly by the user. 
		int n, mm, jsav, cor, dj;
		const double *xx, *yy;
		Base_interp(std::vector<double> &x, const double *y, int m)
			// Constructor: Set up for interpolating on a table of x’s and y’s of length m.Normally called by a derived class, not by the user.
			: n(x.size()), mm(m), jsav(0), cor(0), xx(&x[0]), yy(y)
		{
			dj = std::min(1, (int)std::pow((double)n, 0.25));
		}
		double interp(double x)
		{
			// Given a value x, return an interpolated value, using data pointed to by xx and yy.
			int jlo = cor ? hunt(x) : locate(x);
			return rawinterp(jlo, x);
		}

		int locate(const double x); // See deﬁnitions below.
		int hunt(const double x);

		double virtual rawinterp(int jlo, double x) = 0;
		// Derived classes provide this as the actual interpolation method.

};

struct Linear_interp : Base_interp
	//Piecewise linear interpolation object. Construct with x and y vectors, then call interp for interpolated values. 
{
	Linear_interp(std::vector<double> &xv, std::vector<double> &yv)
		: Base_interp(xv, &yv[0], 2) {}
	double rawinterp(int j, double x) {
		if (xx[j] == xx[j + 1]) return yy[j]; // Table is defective, but we can recover.
		else return yy[j] + ((x - xx[j]) / (xx[j + 1] - xx[j]))*(yy[j + 1] - yy[j]);
	}
};

struct Poly_interp : Base_interp
	// Polynomial interpolation object. Construct with x and y vectors, and the number M of points to be used locally (polynomial order plus one), then call interp for interpolated values.
{
	double dy;
	Poly_interp(std::vector<double> &xv, std::vector<double> &yv, int m)
		: Base_interp(xv, &yv[0], m), dy(0.) {}
	double rawinterp(int jl, double x);
};