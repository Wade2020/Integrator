// intensity integrator.cpp : Defines the entry point for the console application.
// Integrates intensity (cross section) equation from Aharonian into table as well as outputting to files and printing

#include "stdafx.h"
#include <fstream>
#include <random>
#include <math.h>
#include <iostream>
#include <cstdlib>
#include <stdio.h>
#include <iomanip>
#include <ctime>
#include <vector>
#include <algorithm>
#include "Cutpoint.h"
#include "Base_interp.h"
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_integration.h>

using namespace std;

FILE _iob[] = { *stdin, *stdout, *stderr };
extern "C" FILE * __cdecl __iob_func(void) { return _iob; }

int bcounto2;							// Counter for x array (omega 2) for integration
const int nxele = 31;					// Size of x array - number of divisions plus one
vector<double> xvector;					// Defining x array to hold omega 2 values


double drand48()						// Randomly generates a number for use in Monte Carlo between 0 and 1
{
	return rand() / (RAND_MAX + 1.);
}

struct varclass							// Struct to hold variables - all public
{
	long double pi = 3.141592653589793238;
	long double r = 2.8179e-13;				// "r" in equations, in cm.
	double omega1 = 0.01;					// For x ray
	double omega2; //= xvector[bcounto2];		// Setting x values, for gamma ray
	double wavee = omega2;					// "E" in the Aharonian equations
	double mcsqr = 1;						// mc^2 in equations
};

double intensity(double epsilon, long double PI, double rad, double o1, double o2, double E) {		// Cross section (or intensity) equation

	double res =
		(PI*rad*rad / (4 * o1*o1 * o2*o2*o2)) *
		(
		(4 * E*E) / ((E - epsilon)*epsilon) * log((4 * o1 * (E - epsilon)*epsilon) / E)
			- 8 * o1*E
			+ ((2 * (2 * o1*E - 1) * E*E) / ((E - epsilon)*epsilon))
			- ((1 - (1 / (o1*E)))
				* (E*E*E*E / ((E - epsilon) * (E - epsilon) * epsilon*epsilon)))
		);

	return res;
}

double integrand(double x, void *p_param)			// Points to variables and calls eqaution
{
	varclass* p_p = (varclass *)p_param;

	long double P = p_p->pi;
	double ra = p_p->r;
	double om1 = p_p->omega1;
	double om2 = p_p->omega2;
	double En = p_p->wavee;

	return intensity(x, P, ra, om1, om2, En);
}

double int_cs(double om, double Eg, double mc2, double r, long double pi) {		// Calculates the simplified integrated cross section provided by Aharonian
	double z = sqrt(om*Eg) / mc2;			// The omega might be omega 2, or the gamma ray

	double res =

		4 * pi*r*r / pow(z, 4) *
		(
		(z*z - 1 + 1 / (2 * z*z) + log(2 * z)) * log(z + sqrt(z*z - 1))
			+ .5 * log(z)*log(z)
			- .5 * pow(log(z + sqrt(z*z - 1)), 2)
			+ log(2)*log(z)
			- z * sqrt(z*z - 1)
			);

	cout << "Approximation: " << res << "\n";
	return res;
}

int main()
{
	const int nyele = nxele;		// Size of y array - number of divisions plus one
	// xvector corresponds to omega 2 in Aharonian and is defined at the top
	// yarray corresponds to epsilon in Aharonian and is defined below
	// yvector is the same as yarray, but in vector form so the cutpoint method works
	double yarray[nyele];		// Epsilon in Aharonian, corresponding to the energy of the particle
	vector<double> yvector;		// Epsilon in Aharonian, corresponding to the energy of the particle

	//double table[nyele][nxele]; // Table to hold integration values, where order is [row][column]
	vector< vector<double> > table(nxele, vector<double>(nyele));

	double mino2 = 1000;			// X range comes from w1*w2 goes from 10 to 1000 with w1 (omega 1) fixed
	double maxo2 = 100000;
	double stepo2 = log10(maxo2 / mino2) / (nxele - 1);	// Logarithmic step size for omega 2
	double minep, maxep, stepep;
	double abs_error = 1.0e-8;		// Integration quantities required for GSL funciton
	double rel_error = 1.0e-8;
	double abserr = 1;
	size_t neval = 1;
	//void *pointer;
	varclass object;
	gsl_function my_function;		// Defining GSL function to integrate
	my_function.function = &integrand;
	my_function.params = &object;

	int x_if = 10;					// Variable that determines which index of w2 the second analysis corresponds to
									// Used for method of plotting two omega 2 values at once

	bcounto2 = 0;					// Counter for x integration
	while (bcounto2 < nxele)		// Makes the x array, where the values go from 100 to 100000
	{
		xvector.push_back(mino2 * pow(10, bcounto2*stepo2));	// Assinging values to xvector
		object.omega2 = xvector[bcounto2];						// Reassgin value to array
		object.wavee = object.omega2;							// Reassign value to wave ebergy
		//cout << "New omega 2 value: " << xvector[bcounto2] << endl;
		minep = ((object.wavee) / 2)*(1 - sqrt(1 - (1 / (object.omega1 * object.wavee))));		// Sets min and max values for epsilon
		maxep = ((object.wavee) / 2)*(1 + sqrt(1 - (1 / (object.omega1 * object.wavee))));
		//stepep = log10(maxep / minep) / (nyele-1);
		stepep = (maxep - minep) / (nyele - 1);					// Linear step size for epsilon
		//stepep = log10(10 + (pow(10, maxep) - pow(10, minep))) / (nyele-1);
		//pointer = &xvector[bcounto2];
		int bcountep = 0;

		while (bcountep < nyele)	// Makes y array for each value of omega 2 as the max value of y depends on x
		{
			//yarray[bcountep] = maxep / pow(10,stepep * bcountep);
			yarray[bcountep] = maxep - bcountep * stepep;	// Assigns value to y array
			//yarray[bcountep] = maxep / log10(10 + stepep*bcountep);

			// Do the integration now, integrate from 1 to epsilon, save to table value

			gsl_integration_qng(&my_function, minep, yarray[bcountep], abs_error, rel_error, &table[bcounto2][bcountep], &abserr, &neval);

			//cout << "Y array " << yarray[bcountep] << endl;
			bcountep++;
		}

		/*
		cout << "Array of epsilon values for latest omega" << endl;
		for (size_t i = nyele; i--;)
		{
			cout << yarray[i] << endl;
		}
		*/

		
		double index_t_0 = table[bcounto2][0];
		double normal_t_divisor = table[bcounto2][nyele-1] - index_t_0;
		for (int index_count = 0; index_count < nyele; index_count++)
			table[bcounto2][index_count] = (table[bcounto2][index_count] - index_t_0) / normal_t_divisor;

		bcounto2++;
	}
	/*
	cout << "Now printing x array" << endl;
	for (size_t i = nxele; i--;)
	{
		cout << xvector[i] << endl;
	}
	 
	cout << stepo2 <<endl;
	*/

	for (int index = 0; index < nyele; index++) yvector.push_back(yarray[index]);

	double index_x_0 = xvector[0];
	double normal_x_divisor = xvector[nxele - 1] - index_x_0;
	double index_y_0 = yvector[0];
	double normal_y_divisor = yvector[nyele - 1] - index_y_0;
	for (int index_count = 0; index_count < nyele; index_count++) {						// Current loop works only if nxele = nyele; otherwise just use two loops
		xvector[index_count] = (xvector[index_count] - index_x_0) / normal_x_divisor;
		yvector[index_count] = (yvector[index_count] - index_y_0) / normal_y_divisor;
	}

	/*******************************************************************************************************************************************************************************************************/

/*	// For making Monte Carlo simulation
	// Idea is to use for loop to determine if value is between a certain interval in a cutpoint method kind of manner
	// Uses most recent value for omega 2
	// This section is largely unneeded due to the cutpoint code from Dr. Timokhin

	double rand_num;					// The random number itself
	const int nrand_num = 10000;		// Number of random numbers
	vector<double> rand_val;			// Array to hold the value corresponding to the random number in the integration table
	srand(time(0));						// Attempt to "randomly" seed values with seconds value
	for (int rand_count = 0; rand_count < nrand_num; rand_count++)	// Test several random numbers
	{
		rand_num = drand48();

		for (int rand_tab_count = 1; rand_tab_count < nyele; rand_tab_count++)	// Loops to find where the random number is	
			if (rand_num < table[nxele - 1][rand_tab_count])	// If the random value is less than a value of the table...
			{
				rand_val.push_back(yvector[rand_tab_count]);			// ... than that value of epsilon is recorded for later histogram plotting...
				break;													// ... and the code continues on; the value corresponding to the random number is found
			}												// Otherwise the loop continues trying to find the value
	}

	// Plots the distribution of integrated values to get a sense of the probability outcome
	//for (int rand_tab_count = 0; rand_tab_count < nyele; rand_tab_count++)
		//cout << table[nxele - 1][rand_tab_count] << endl;
	
	vector<double> rand_val2;			// Array to hold the value coressponding to the random number in the integration table
	for (int rand_count = 0; rand_count < nrand_num; rand_count++)	// Test several random numbers
	{
		rand_num = drand48();

		for (int rand_tab_count = 1; rand_tab_count < nyele; rand_tab_count++)	// Loops to find where the random number is
			if (rand_num < table[x_if][rand_tab_count])	// If the random value is less than a value of the table...
			{
				rand_val2.push_back(yvector[rand_tab_count]);			// ... than that value of epsilon is recorded for later histogram plotting...
				break;													// ... and the code continues on; the value corresponding to the random number is found
			}												// Otherwise the loop continues trying to find the value
	}
	// Plots the distribution of integrated values to get a sense of the probability outcome
	//for (int rand_tab_count = 0; rand_tab_count < nyele; rand_tab_count++)
		//cout << table[x_if][rand_tab_count] << endl;

	//cout << "stepep " << stepep << " Max " << maxep << " Min " << minep << endl;	// Outputs epsilon values for reference in Python plotting

	// Output some values, here a set of values corresponding to Monte Carlo methods, perhaps using cutpoint for one value of omega 2 (the last one), for plotting in Python
/*	outputFile.open("MC.txt");			// Output to text file called "MC"

	for (int output_counter = 0; output_counter < rand_num_num - 1; output_counter++) {
		outputFile << rand_value[output_counter] << endl;		// Output counter value (not needed) and value assigned by randomization
	}
	outputFile << rand_value[rand_num_num - 1];				// Does the same as above, just without endlining eliminate index error in Python

	outputFile.close();					// Stop outputting to file

	
	/*******************************************************************************************************************************************************************************************************/
	
	// Uses cutpoint for Monte Carlo simulations 

	Cutpoint cm;
	vector<double> table_column(nyele);							// Array for different epsilons for interpolated omega 2 (last one)
	vector<double> table_column2(nyele);						// Array for different epsilons for interpolated omega 2 (x_if one)
	int ncp = (nyele - 1) / 2;									// Number of cutpoints
			
	// Saves the most recent table column (for the last set of epsilons and last omega 2 value)
	for (int t_c_c = 0; t_c_c < nyele; t_c_c++)
	{
		table_column[t_c_c] = table[nxele - 1][t_c_c];
		table_column2[t_c_c] = table[x_if][t_c_c];
	}

	table_column[0] = 0.001;
	cm.Initialize(yvector, table_column, ncp);
	//const int nksi = 6;
	//double ksi[nksi] = { .72, .36, .94, .37, 0.0002, 0 };
	//double cmksi[nksi];
	const int nrand_num = 10000;				// Number of random numbers
	vector<double> ksi(nrand_num), rand_val(nrand_num), rand_val2(nrand_num);		// Array of random numbers and corresponding cutpoints
	srand(time(0));						// Attempt to "randomly" seed values with seconds value

	for (int i = 0; i<nrand_num; i++)
	{
		ksi[i] = drand48();
		rand_val[i] = cm(ksi[i]);
		//cout << "ksi=" << ksi[i] << "  x=" << rand_val[i] << endl;
	}
	
	cm.Initialize(yvector, table_column2, ncp);
	table_column2[0] = 0.001;

	for (int i = 0; i<nrand_num; i++)
	{
		rand_val2[i] = cm(ksi[i]);
		//cout << "ksi=" << ksi[i] << "  x=" << rand_val2[i] << endl;
	}


	/*******************************************************************************************************************************************************************************************************/

	// Make table of cutpoint indices, find value of omega 2 in the xvector, then cutpoint and output

	// Table of indices and array to hold column
	vector< vector<int> > cp_table(nxele, vector<int>(ncp));
	vector<int> tab_col(ncp);
	vector<int> effi_vec;

	for (int i = 0; i < nxele; i++) {
		for (int t_c_c = 0; t_c_c < nyele; t_c_c++) {
			table_column[t_c_c] = table[i][t_c_c];
		};

		table_column[0] = 0.001;
		cm.Initialize(yvector, table_column, ncp);

		for (int j = 0; j < ncp; j++)
			cp_table[i][j] = cm.return_I(j);
	}

	// Pick the omega value and find it
	double om2_val = 0.5;							// Some value for normalized omega 2
	cm.Initialize(xvector, xvector, (nxele-1) / 2);
	cm(om2_val);									// Finds closest value of omega 2 to value put in
	int om2_ind = cm.return_z();					// om2ind could be replaced by or replace x_if


	// Now use cutpoint for the desired value of omega 2
	for (int t_c_c = 0; t_c_c < nyele; t_c_c++)
		table_column[t_c_c] = table[om2_ind][t_c_c];

	for (int t_c_c = 0; t_c_c < ncp; t_c_c++)
		tab_col[t_c_c] = cp_table[om2_ind][t_c_c];

	for (int i = 0; i<nrand_num; i++)
	{
		rand_val[i] = cm(ksi[i], tab_col, yvector, table_column, ncp);
		effi_vec.push_back(cm.co);			// Part of efficiency calculation
		//cout << i << endl;
		//cout << "ksi=" << ksi[i] << "  x=" << rand_val[i] << " Efficiency: " << cm.co << endl;
	}

	// Compute average efficiency
	double sum = 0;
	for (int x : effi_vec) sum += x;
	cout << "Average Efficiency: " << sum / effi_vec.size() << endl;
	
	ofstream outputFile;
	outputFile.open("MC.txt");						// Output to text file called "MC"
	
	for (int output_counter = 0; output_counter < nyele; output_counter++) {
		outputFile << yvector[output_counter] << " " << table[om2_ind][output_counter] << " " << rand_val[output_counter] << endl;	// Output y array value and integration value
	}

	for (int output_counter = nyele; output_counter < (nrand_num - 1); output_counter++) {
		outputFile << 0 << " " << 0 << " " << rand_val[output_counter] << endl;	// Output y array value and integration value
	}
	outputFile << 0 << " " << 0 << " " << rand_val[nrand_num - 1];	// Output y array value and integration value
	outputFile.close();					// Stop outputting to file


	/*******************************************************************************************************************************************************************************************************/

	// Output some values, here just one column of the integration table for one value of omega 2 (the last one), for plotting in Python
	// Current loop only works if nxele = nyele

/*	ofstream outputFile;
	outputFile.open("Data.txt");		// Output to text file called "Data"
	
	for (int output_counter = 0; output_counter < nyele; output_counter++) {
	outputFile << yvector[output_counter] << " " << table[nxele - 1][output_counter] << " " << rand_val[output_counter] << " " << table[x_if][output_counter] << " " << rand_val2[output_counter] << endl;	// Output y array value and integration value
	}

	for (int output_counter = nyele; output_counter < (nrand_num-1); output_counter++) {
		outputFile << 0 << " " << 0 << " " << rand_val[output_counter] <<  " " << 0 << " " << rand_val2[output_counter] << endl;	// Output y array value and integration value
	}
	outputFile << 0 << " " << 0 << " " << rand_val[nrand_num-1] << " " << 0 << " " << rand_val2[nrand_num-1];	// Output y array value and integration value
	outputFile.close();					// Stop outputting to file
	
	/*for (int output_counter = 0; output_counter < nyele; output_counter++) {
		outputFile << yvector[output_counter] << " " << table[nxele - 1][output_counter] << endl;	// Output y array value and integration value
	}
	outputFile << yvector[nyele - 1] << " " << table[nxele - 1][nyele - 1];				// Does the same as above, just without endlining eliminate index error in Python
	outputFile.close();					// Stop outputting to file

	/*
	//Outputs values of cross section equation (function is correct); if you use, remember to close file
	for (double output_counter = minep; output_counter < maxep; output_counter++) {
	outputFile << output_counter << " " << intensity(output_counter, object.pi, object.r, object.omega1, object.omega2, object.wavee) << endl;
	}
	

	/*******************************************************************************************************************************************************************************************************/

	// Outputs table values

	// Output format taken from Ohio State Example (https://www.physics.ohio-state.edu/~ntg/780/gsl_examples/qags_test.cpp)

	//cout.setf(ios::fixed, ios::floatfield);	// output in fixed format
	/*cout.precision(5);		// 5 digits in doubles

	int width = 15;  // setw width for output
	bool print_val = false;
	int mod_step = nyele / 6;						// Set up for nxele = nyele
												// For only outputting so many values in the table
	
	cout << "Now outputing the table values with every " << mod_step << " shown" << endl;

	/*
	cout << table[0][0] << setw(width) << table[0][1] << setw(width) << table[0][2] << setw(width) << table[0][3] << setw(width) << table[0][4] << endl;
	cout << table[1][0] << setw(width) << table[1][1] << setw(width) << table[1][2] << setw(width) << table[1][3] << setw(width) << table[1][4] << endl;
	cout << table[2][0] << setw(width) << table[2][1] << setw(width) << table[2][2] << setw(width) << table[2][3] << setw(width) << table[2][4] << endl;
	
	cout << endl << endl << endl;
	*/	

	/*for (int y_print_count = 0; y_print_count < nyele; y_print_count++) {					// Loop over rows of table
		for (int x_print_count = 0; x_print_count < nxele; x_print_count++) {				// But first loop over columns
			if (((y_print_count % mod_step) == 1) && ((x_print_count % mod_step) == 1))		// Print every mod_step time
			{
				cout << table[x_print_count][y_print_count] << setw(width);					// Print integration values
				print_val = true;															// Document that the print occured
			}
		}
		if (print_val) cout << table[nxele-1][y_print_count] << endl;					// Print without space to align next line correctly only 
		print_val = false;
	}

	cout << endl;

	// Test more carefully
	for (size_t j = nyele; j--;){
		for (size_t i = nxele; i--;){ 
			if (((j % mod_step) == 1) && ((i % mod_step) == 1))		// Print every mod_step time
			{
				cout << table[i][j] << setw(width);					// Print integration values
				print_val = true;															// Document that the print occured
			}
		}
		if (print_val) cout << table[nxele - 1][j] << endl;					// Print without space to align next line correctly only 
		print_val = false;
	}	
	
	
	cout << table[0][0] << setw(width) << table[1][0] << setw(width) << table[2][0] << setw(width) << table[3][0] << setw(width) << table[4][0] << endl;
	cout << table[0][1] << setw(width) << table[1][1] << setw(width) << table[2][1] << setw(width) << table[3][1] << setw(width) << table[4][1] << endl;
	cout << table[0][2] << setw(width) << table[1][2] << setw(width) << table[2][2] << setw(width) << table[3][2] << setw(width) << table[4][2] << endl;


	/*******************************************************************************************************************************************************************************************************/
	
/*	// Interpolates values for the table, per epsilon, based on a value given for omega 2. Uses the Numerical Recipes code at the top of the program.

	double rand_num;					// The random number itself
	//const int nrand_num = 10000;		// Number of random numbers
	//vector<double> table_column;							// Array for different epsilons for interpolated omega 2
	double om2_interp = 1800;								// Value of omega 2 to interpolate for
	cout << endl << "Now outputting interpolation values for omega 2 equals " << om2_interp << endl;

	for (int interp_y_count = 0; interp_y_count < nyele; interp_y_count++)
	{
		vector<double> table_row;			// Holding the row of the table being interpolated
											// The "xx" vector is xvector, "yy" array is table_row

		for (int interp_x_count = 0; interp_x_count < nxele; interp_x_count++)
		{
			table_row.push_back(table[interp_x_count][interp_y_count]);		// Reassinging each row of the table to an array to interpolate from
		};

		// For linear, use Linear_interp myfunc(xvector, table_row), for polynomial use Poly_interp myfunc(xvector,table_row,4), where 4 is cubic (or just do any #) 
		Linear_interp myfunc(xvector, table_row);
		table_column.push_back(myfunc.interp(om2_interp));		// Makes the array of table values for each epsilon
	}

	// Now Monte Carlo and plot

	vector<double> rand_value3;			// Array to hold the value coressponding to the random number in the integration table
	for (int rand_count = 0; rand_count < nrand_num; rand_count++)	// Test several random numbers
	{
		rand_num = drand48();

		// Normailizes if not already done so
		// Maybe put in function to normalize stepwise (so I don't have same code a few times)
		for (int rand_tab_count = 1; rand_tab_count < nyele; rand_tab_count++)	// Loops to find where the random number is
		{
			if (rand_num < table_column[rand_tab_count])	// If the random value is less than a value of epsilon...
			{
				rand_value3.push_back(yvector[rand_tab_count]);		// ... than that value of epsilon is recorded for later histogram plotting...
				break;													// ... and the code continues on; the value corresponding to the random number is found
			}													// Otherwise the loop continues trying to find the value
		}
	}

/*	// Plots the distribution of integrated values to get a sense of the probability outcome
	// Could put the output statement right after te assignment to reduce lines, but this is a nice part on its own
	cout << endl << "Now outputting interpolated, normalized table values" << endl;
	for (int rand_tab_count = 0; rand_tab_count < nyele; rand_tab_count++) {
		cout << table_column[rand_tab_count] << endl;
	}

	fstream outputFile;		// Delete when finished


	outputFile.open("Interp.txt");		// Output to text file called "Interp"
	for (int output_counter = 0; output_counter < nyele; output_counter++) {
		outputFile << yvector[output_counter] << " " << table[nxele - 1][output_counter] << " " << table_column[output_counter] << endl;	// Output y array value and integration value
	}
	outputFile.close();					// Stop outputting to file

	/*******************************************************************************************************************************************************************************************************/
	/*
	// Test GSL function
	double x, y;
	x = 5.0;
	y = gsl_sf_bessel_J0(x);
	printf("J0(%g) = %.18e\n", x, y);
	*/

	// My way of keeping window open (for Windows machine)
	cout << "Done";
	int num;
	cin >> num;
	return 0;
}