// Parameter estimation for GARCH(1,1) model
// Fitted to index levels of a hypothetical asset C

#include "stdafx.h"
#using <mscorlib.dll>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <algorithm>
using namespace System;
using namespace std;

// Function to calculate the mean value of
// a set of n vectors each of dimension n
// namely a (n x n) matrix
vector<double> VMean(vector<vector<double> > X, int n)
{
    vector<double> meanX(n);
	for (int i=0; i<=n-1; i++)
	{
		meanX[i]=0.0;
		for (int j=0; j<=n-1; j++)
            {
				meanX[i] += X[i][j];
			}
		meanX[i] = meanX[i] / n;
	}
return meanX;
}

// Function to add two vectors together
vector<double> VAdd(vector<double> x, vector<double> y)
{
    int n = x.size();
    vector<double> z(n, 0.0);
	for (int j=0; j<=n-1; j++)
		z[j] = x[j] + y[j];
return z;
}

// Function to subtract two vectors
vector<double> VSub(vector<double> x, vector<double> y)
{
    int n = x.size();
    vector<double> z(n, 0.0);
	for (int j=0; j<=n-1; j++)
		z[j] = x[j] - y[j];
return z;
}

// Function to multiply a vector by a constant
vector<double> VMult(vector<double> x, double a)
{
    int n = x.size();
    vector<double> z(n, 0.0);
	for (int j=0; j<=n-1; j++)
		z[j] = a*x[j];
return z;
}

// Sum of a vector
double VecSum(vector<double> x)
{
	int n = x.size();
	double Sum = 0.0;
	for (int i=0; i<=n-1; i++)
		Sum += x[i];
	return Sum;
}

// Calculates unbiased sample variance
double VecVar(vector<double> x)
{
	double n = x.size();
	double sumM = 0.0;
	for (int i=0; i<=n-1; i++)
		sumM += x[i];
	double mean = sumM / n;
	double sumV = 0.0;
	for (int i=0; i<=n-1; i++)
		sumV += (x[i] - mean)*(x[i] - mean);
	return sumV / (n-1);
}

// Nelder Mead Algorithm
vector<double> NelderMead(double (*f)(vector<double>), int N, double NumIters,double MaxIters,
						  double Tolerance, vector<vector<double> > x)
{
int i,j;

// Value of the function at the vertices
vector<vector<double> > F(N+1, vector<double>(2));

// Step 0.  Ordering and Best and Worst points
// Order according to the functional values, compute the best and worst points
step0:
NumIters = NumIters + 1;
for (j=0; j<=N; j++){
    vector<double> z(N, 0.0);								// Create vector to contain
	for (i=0; i<=N-1; i++)
		z[i] = x[i][j];
	F[j][0] = f(z);		                 					// Function values
	F[j][1] = j;											// Original index positions
}
sort(F.begin(), F.end());

// New vertices order first N best initial vectors and
// last (N+1)st vertice is the worst vector
// y is the matrix of vertices, ordered so that the worst vertice is last
vector<vector<double> > y(N, vector<double>(N+1));
for (j=0; j<=N; j++) {
	for (i=0; i<=N-1; i++) {
		y[i][j] = x[i][F[j][1]];
	}
}

//  First best vector y(1) and function value f1
vector<double> x1(N, 0.0); for (i=0; i<=N-1; i++) x1[i] = y[i][0];
double f1 = f(x1);

// Last best vector y(N) and function value fn
vector<double> xn(N, 0.0); for (i=0; i<=N-1; i++) xn[i] = y[i][N-1];
double fn = f(xn);

// Worst vector y(N+1) and function value fn1
vector<double> xn1(N, 0.0); for (i=0; i<=N-1; i++) xn1[i] = y[i][N];
double fn1 = f(xn1);

// z is the first N vectors from y, excludes the worst y(N+1)
vector<vector<double> > z(N, vector<double>(N));
for (j=0; j<=N-1; j++) {
	for (i=0; i<=N-1; i++)
		z[i][j] = y[i][j];
}

// Mean of best N values and function value fm
vector<double> xm(N, 0.0); xm = VMean(z,N);
double fm = f(xm);

// Reflection point xr and function fr
vector<double> xr(N, 0.0); xr = VSub(VAdd(xm, xm), xn1);
double fr = f(xr);

// Expansion point xe and function fe
vector<double> xe(N, 0.0); xe = VSub(VAdd(xr, xr), xm);
double fe = f(xe);

// Outside contraction point and function foc
vector<double> xoc(N, 0.0);	xoc = VAdd(VMult(xr, 0.5), VMult(xm, 0.5));
double foc = f(xoc);

// Inside contraction point and function foc
vector<double> xic(N, 0.0);	xic = VAdd(VMult(xm, 0.5), VMult(xn1, 0.5));
double fic = f(xic);

while ((NumIters <= MaxIters) && (abs(f1-fn1) >= Tolerance))
{
// Step 1. Reflection Rule
if ((f1<=fr) && (fr<fn)) {
    for (j=0; j<=N-1; j++) {
		for (i=0; i<=N-1; i++)  x[i][j] = y[i][j]; }
		for (i=0; i<=N-1; i++)	x[i][N] = xr[i];
	goto step0;
}

// Step 2.  Expansion Rule
if (fr<f1) {
	for (j=0; j<=N-1; j++) {
		for (i=0; i<=N-1; i++)  x[i][j] = y[i][j]; }
	if (fe<fr)
		for (i=0; i<=N-1; i++)	x[i][N] = xe[i];
	else
		for (i=0; i<=N-1; i++)	x[i][N] = xr[i];
	goto step0;
}

// Step 3.  Outside contraction Rule
if ((fn<=fr) && (fr<fn1) && (foc<=fr)) {
	for (j=0; j<=N-1; j++) {
		for (i=0; i<=N-1; i++)  x[i][j] = y[i][j]; }
		for (i=0; i<=N-1; i++)	x[i][N] = xoc[i];
	goto step0;
}

// Step 4.  Inside contraction Rule
if ((fr>=fn1) && (fic<fn1)) {
	for (j=0; j<=N-1; j++) {
		for (i=0; i<=N-1; i++)  x[i][j] = y[i][j]; }
		for (i=0; i<=N-1; i++)	x[i][N] = xic[i];
	goto step0;
}

// Step 5. Shrink Step
for (i=0; i<=N-1; i++)
	x[i][0] = y[i][0];
for (i=0; i<=N-1; i++) {
    for (j=1; j<=N; j++)
		x[i][j] = 0.5*(y[i][j] + x[i][0]);
}
	goto step0;
}

// Output component
vector<double> out(N+2);
for (i=0; i<=N-1; i++)
	out[i] = x1[i];
	out[N] = f1;
	out[N+1] = NumIters;
return out;

}
// Log Likelihood for GARCH(1,1)
// Returns the negative Log Likelihood, for minimization
// Variance equation is
// h(t) = B[0] + B[1]*ret2[i+1] + B[2]*h[i+1];
// where B[0]=omega, B[1]=alpha, B[2]=beta

double LogLikelihood(vector<double> B)
{
	ifstream inPrices;
	inPrices.open("AssetC_IndexPrices.txt");					// Text file with AssetC_Index Price levels
	vector<double> Price(0, 0.0);
	int i=0;
	double n1;
	while (!inPrices.eof()) {
		inPrices >> n1;
		Price.push_back(n1);
	i++ ;
	}
	inPrices.close();
	int n = Price.size();
	vector<double> ret(n-1);					// Return
	vector<double> ret2(n-1);					// Return squared
	vector<double> GARCH(n-1, 0.0);				// GARCH(1,1) variance
	vector<double> LogLike(n-1, 0.0);
												// Penalty for non-permissible parameter values
	if  ( (B[0]<0.0) || (B[1]<0.0) || (B[2]<0.0) ||(B[1]+B[2]>=1) )
		return 1e100;
	else										// Construct the log likelihood
	for (i=0; i<=n-2; i++) {
		ret[i] = log(Price[i]/Price[i+1]);
		ret2[i] = ret[i]*ret[i];
	}
	GARCH[n-2] = VecVar(ret);
	LogLike[n-2] = -log(GARCH[n-2]) - ret2[n-2]/GARCH[n-2];
	for (i=n-3; i>=0; i--) {
		GARCH[i] = B[0] + B[1]*ret2[i+1] + B[2]*GARCH[i+1];
	  LogLike[i] = -log(GARCH[i]) - ret2[i]/GARCH[i];
	}
	return -VecSum(LogLike);
}

int main()
{
	int N = 3;						// Number of GARCH(1,1) parameters
	double NumIters = 1;			// First Iteration
	double MaxIters = 1e2 ;			// Maximum number of iterations
	double Tolerance = 1e-15;		// Tolerance on best and worst function values

	vector<vector<double> > s(N, vector<double>(N+1));
	// Vertice 0	Vertice 1		Vertice 2		Vertice 3		Vertice 4
	s[0][0]=0.00002;	s[0][1]=0.00001;	s[0][2]=0.00000015;		s[0][3]=0.0000005;
	s[1][0]=0.10;		s[1][1]=0.11;		s[1][2]=0.09;			s[1][3]=0.15;
	s[2][0]=0.85;		s[2][1]=0.87;		s[2][2]=0.90;			s[2][3]=0.83;
	vector<double> NM = NelderMead(LogLikelihood, N, NumIters, MaxIters, Tolerance,s);
	cout << "GARCH(1,1) Parameters" << endl;
	cout << "--------------------------------" << endl;
	cout << "Omega       = " << NM[0] << endl;
	cout << "Alpha       = " << NM[1] << endl;
	cout << "Beta        = " << NM[2] << endl;
	cout << "Persistence = " << NM[1] + NM[2] << endl;
	cout << "--------------------------------" << endl;
	cout << "Log Likelihood value = " << NM[3] << endl;
	cout << "Number of iterations = " << NM[4] << endl;
	return 0;
}
