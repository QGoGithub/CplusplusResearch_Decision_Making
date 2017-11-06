// Black Scholes implied volatility using S&P500 option prices
#include "stdafx.h"
#using <mscorlib.dll>
#include <complex>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <algorithm>
#include <time.h>
using namespace System;
using namespace std;

// Cumulative Standard Normal
double N(double x) {
	const double pi =  4.0*atan(1.0);
    const double b1 =  0.319381530;
    const double b2 = -0.356563782;
	const double b3 =  1.781477937;
	const double b4 = -1.821255978;
	const double b5 =  1.330274429;
	double t = 1.0 / (1.0 + 0.2316419*x);
	double n = exp(-x*x/2.0)/sqrt(2*pi);
	return 1.0 - n*(b1*t + b2*pow(t,2) + b3*pow(t,3) + b4*pow(t,4) + b5*pow(t,5));
}
// Black-Scholes Price
double BSPrice(double S, double K, double r, double T, double v, char PutCall) {
	double d = (log(S/K) + T*(r + 0.5*v*v)) / (v*sqrt(T));
	double Call = S*N(d) - exp(-r*T)*K*N(d - v*sqrt(T));
	if (PutCall=='C')
		return Call;
	else
	return Call - S + K*exp(-r * T);
}
// Bisection Algorithm
double BisecBSV(double S, double K, double r, double T,
				double a, double b, double MktPrice, char PutCall) {
	const int MaxIter = 500;
	double Tol = 0.00001;
	double midP, midCdif;
	double  lowCdif = MktPrice - BSPrice(S, K, r, T, a, PutCall);
	double highCdif = MktPrice - BSPrice(S, K, r, T, b, PutCall);
	if (lowCdif*highCdif > 0)
		return -1;
	else
	for (int i=0; i<=MaxIter; i++) {
	    midP = (a + b) / 2.0;
		midCdif = MktPrice - BSPrice(S, K, r, T, midP, PutCall);
		if (abs(midCdif)<Tol) goto LastLine;
			else {
		        if (midCdif>0) a = midP;
			else b = midP;
			}
		}
	LastLine:
	return midP;
}

// Main program
int main() {
	double S = 91.04;				// Spot Price
	double T = 0.1205;				// Maturity in Years
	double r = 0.0025;				// Interest Rate
	double a = 0.000001;			// Bisection algorithm starting value
	double b = 10.0;				// Bisection algorithm starting value

	// Open the option prices for SPY
	ifstream inPrices;
	inPrices.open("SPY.txt");
	vector<double> PutCall(0, 0.0);			// 1=Call 2=Put
	vector<double> K(0, 0.0);				// Strike
	vector<double> MPrice(0, 0.0);			// Market Price
	double n1, n2, n3;
	int i=0;
	while (!inPrices.eof()) {
		inPrices >> n1 >> n2 >> n3;
		PutCall.push_back(n1);
		K.push_back(n2);
		MPrice.push_back(n3);
	i++ ;
	}
	int n = PutCall.size();
	vector<double> IV(n);
	cout << "Put/Call    Strike    Price     B-S Implied Vol " << endl;
	for (int i=0; i<n; i++) {
		IV[i] = BisecBSV(S, K[i], r, T, a, b, MPrice[i], PutCall[i]);
		cout << PutCall[i] << "              " << K[i] << "      " << MPrice[i] << "      " << IV[i] << endl;
	}
	return 0;
}
