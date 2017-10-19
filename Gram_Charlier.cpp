Option price - Gram Charlier


#include "stdafx.h"
#using <mscorlib.dll>
using namespace System;
using namespace std;
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>

// Standard Normal density f(x)
double fz(double x) {
	double pi = 4.0*atan(1.0);
	return exp(-x*x*0.5)/sqrt(2*pi);
}

// Boole's Rule
double Boole(double StartPoint, double EndPoint, int n) {
vector<double> X(n+1, 0.0);
vector<double> Y(n+1, 0.0);
double delta_x = (EndPoint - StartPoint)/double(n);
for (int i=0; i<=n; i++) {
	X[i] = StartPoint + i*delta_x;
	Y[i] = fz(X[i]);
}
double sum = 0;
for (int t=0; t<=(n-1)/4; t++) {
	int ind = 4*t;
    sum += (1/45.0)*(14*Y[ind] + 64*Y[ind+1] + 24*Y[ind+2] + 64*Y[ind+3] + 14*Y[ind+4])*delta_x;
}
return sum;
}
// Cumulative Standard Normal distribution F(x)
double N(double x) {
	return Boole(-10, x, 240);
}

// Gram Charlier option price
double GramCharlier(double S, double K, double r, double q, double T,
					double v, double skew, double kurt, char PutCall) {
	double r1 = r/12;
	double v1 = v/sqrt(12.0);
	double Nskew = skew / sqrt(T);
	double Nkurt = kurt / T;
	double Nvol = sqrt(T)*v1;
	double d = (log(S/K) + T*(r1 - q) + Nvol*Nvol/2) / Nvol;
	double Call = S*exp(-q*T)*N(d) - K*exp(-r1*T)*N(d - Nvol) +
	              S*exp(-q*T)*fz(d)*Nvol*
	             ((Nskew/6)*(2*Nvol - d) - (Nkurt/24)*(1 - d*d + 3*d*Nvol - 3*Nvol*Nvol));
	if (PutCall=='C')
		return Call;
	else
		return Call + K*exp(-r1*T) - S*exp(-q*T);
}

int main() {
	double Spot = 30.0;			// Spot of the underlying
	double Strike = 30.0;		// Strike price
	double r = 0.05;			// Risk free rate
	double q = 0.0;				// Foreign Risk free Rate
	double T = 5;				// Years to maturity
	double v = 0.30;			// Volatility of the underlying
	double skew = -2.3;			// Skewness of the underlying
	double kurt = 1.2;			// Kurtosis of the underlying
	cout << "The Gram-Charlier Call Price is " << GramCharlier(Spot, Strike, r, q, T, v, skew, kurt, 'C') << endl;
	cout << "The Gram-Charlier Put Price  is " << GramCharlier(Spot, Strike, r, q, T, v, skew, kurt, 'P') << endl;
}
