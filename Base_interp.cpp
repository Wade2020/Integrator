#include "stdafx.h"
#include <vector>
#include <stdio.h>
#include <algorithm>
#include "Base_interp.h"

using namespace std;

// Pulling this from Numerical Recipes

int Base_interp::locate(const double x)
// Given a value x, return a value j such that x is(insofar as possible) centered in the subrange xx[j..j + mm - 1], where xx is the stored pointer.The values in xx must be monotonic, either increasing or decreasing.The returned value is not less than 0, nor greater than n - 1.
{
	int ju, jm, jl;
	if (n < 2 || mm < 2 || mm > n)throw("locate size error");
	bool ascnd = (xx[n - 1] >= xx[0]);				// True if ascending order of table, false otherwise.
	jl = 0;											// Initialize lower 
	ju = n - 1;										// and upper limits.
	while (ju - jl > 1) {							// If we are not yet done, 
		jm = (ju + jl) >> 1;						// compute a midpoint, 
		if (x >= xx[jm] == ascnd)
			jl = jm;								// and replace either the lower limit
		else ju = jm;								// or the upper limit, as appropriate. 
	}												// Repeat until the test condition is satisﬁed.
	cor = abs(jl - jsav) > dj ? 0 : 1;				// Decide whether to use hunt or locate next time.
	jsav = jl;
	return max(0, min(n - mm, jl - ((mm - 2) >> 1)));
}

int Base_interp::hunt(const double x) // Given a value x, return a value j such that x is(insofar as possible) centered in the subrange xx[j..j + mm - 1], where xx is the stored pointer. The values in xx must be monotonic, either increasing or decreasing. The returned value is not less than 0, nor greater than n - 1. 
{
	int jl = jsav, jm, ju, inc = 1;
	if (n < 2 || mm < 2 || mm > n) throw("hunt size error");
	bool ascnd = (xx[n - 1] >= xx[0]);						// True if ascending order of table, false otherwise. 
	if (jl < 0 || jl >n - 1) {								// Input guess not useful.Go immediately to bisection.
		jl = 0; ju = n - 1;
	}
	else {
		if (x >= xx[jl] == ascnd) {							// Hunt up: 
			for (;;) {
				ju = jl + inc;
				if (ju >= n - 1) { ju = n - 1; break; }		// Oﬀ end of table. 
				else if (x < xx[ju] == ascnd) break;		// Found bracket. 
				else {										// Not done, so double the increment and try again.
					jl = ju;
					inc += inc;
				}
			}
		}
		else {												// Hunt down : 
			ju = jl;
			for (;;) {
				jl = jl - inc;
				if (jl <= 0) { jl = 0; break; }				// Oﬀ end of table. 
				else if (x >= xx[jl] == ascnd) break;		// Found bracket. 
				else {										// Not done, so double the increment and try again.
					ju = jl;
					inc += inc;
				}
			}
		}
	}

	while (ju - jl > 1) {				// Hunt is done, so begin the ﬁnal bisection phase : 
		jm = (ju + jl) >> 1;
		if (x >= xx[jm] == ascnd)
			jl = jm;
		else ju = jm;
	}
	cor = abs(jl - jsav) > dj ? 0 : 1;	// Decide whether to use hunt or locate next time.
	jsav = jl;
	return min(0, min(n - mm, jl - ((mm - 2) >> 1)));
}

double Poly_interp::rawinterp(int jl, double x)
// Given a value x, and using pointers to data xx and yy, this routine returns an interpolated value y, and stores an error estimate dy. The returned value is obtained by mm - point polynomial interpolation on the subrange xx[jl..jl + mm - 1]. 
{
	int i, m, ns = 0;
	double y, den, dif, dift, ho, hp, w;
	const double *xa = &xx[jl], *ya = &yy[jl];
	vector<double> c(mm), d(mm);
	dif = abs(x - xa[0]);
	for (i = 0; i < mm; i++) { // Here we ﬁnd the index ns of the closest table entry,
		if ((dift = abs(x - xa[i])) < dif) {
			ns = i;
			dif = dift;
		}
		c[i] = ya[i]; // and initialize the tableau of c’s and d’s.
		d[i] = ya[i];
	}
	y = ya[ns--]; // This is the initial approximation to y. 
	for (m = 1; m < mm; m++) { //For each column of the tableau,
		for (i = 0; i < mm - m; i++) { // we loop over the current c’s and d’s and update them.
			ho = xa[i] - x;
			hp = xa[i + m] - x;
			w = c[i + 1] - d[i];
			if ((den = ho - hp) == 0.0) throw("Poly_interp error"); // This error can occur only if two input xa’s are(to within roundoﬀ) identical.
			den = w / den;
			d[i] = hp*den; // Here the c’s and d’s are updated.
			c[i] = ho*den;
		}
		y += (dy = (2 * (ns + 1) < (mm - m) ? c[ns + 1] : d[ns--]));
		// After each column in the tableau is completed, we decide which correction, c or d, we want to add to our accumulating value of y, i.e., 
		// which path to take through the tableau — forking up or down. We do this in such a way as to take the most “straight line” route 
		// through the tableau to its apex, updating ns accordingly to keep track of where we are. This route keeps the partial approximations
		// centered (insofar as possible) on the target x.The last dy added is thus the error indication.
	}
	return y;
}