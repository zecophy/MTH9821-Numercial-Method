#pragma once
#include"Generator.h"

class BasketOptionPricer {
public:
	int N;
	Option op1;
	Option op2;
	double rho;
	BoxMuller gen;

	vector<double> Z1;
	vector<double> Z2;
	vector<double> S1;
	vector<double> S2;
	vector<double> Vs;

	double value;

	BasketOptionPricer(Option op1_,Option op2_,double rho_)
		:op1(op1_),op2(op2_),rho(rho_){}

	void generateZs() {
		int n = Z1.size();
		while (n < N) {
			Z1.push_back(gen());
			Z2.push_back(gen());
			n += 1;
		}
		return;
	}

	void generateSs() {
		generateZs();
		int n = S1.size();
		while (n < N) {
			double s1 = op1.S0*exp((op1.r - op1.q - 0.5*op1.sigma*op1.sigma)*op1.T 
				+ op1.sigma*sqrt(op1.T)*Z1[n]);
			double s2 = op2.S0*exp((op2.r - op2.q - 0.5*op2.sigma*op2.sigma)*op2.T
				+ op2.sigma*sqrt(op2.T)*(rho*Z1[n] + sqrt(1 - rho * rho)*Z2[n]));
			S1.push_back(s1);
			S2.push_back(s2);
			n += 1;
		}
		return;
	}

	void generateVs() {
		generateSs();
		int n = Vs.size();
		while (n < N) {
			double v = exp(-op1.r*op1.T)*max(S1[n] + S2[n] - op1.K, 0.0);
			Vs.push_back(v);
			n += 1;
		}
		return;
	}

	void PricerGo(int N_) {
		N = N_;
		generateVs();
		value = accumulate(Vs.begin(), Vs.end(), 0.0) / (double)N;
		return;
	}



};
