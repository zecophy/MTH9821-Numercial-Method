#pragma once
#include"Generator.h"

class BasketOptionPathDependentPricer {
public:
	int N, M;
	Option op1;
	Option op2;
	double rho;
	BoxMuller gen;

	vector<vector<double>> Z1;
	vector<vector<double>> Z2;
	vector<vector<double>> S1;
	vector<vector<double>> S2;
	vector<double> Vs;

	double value;

	BasketOptionPathDependentPricer(Option op1_, Option op2_, double rho_, double m)
		:op1(op1_), op2(op2_), rho(rho_), M(m) {}

	void generateZs() {
		int n = Z1.size();
		while (n < N) {
			vector<double> temp1;
			vector<double> temp2;
			for (int i = 0;i < M;i++) {
				temp1.push_back(gen());
				temp2.push_back(gen());
			}
			Z1.push_back(temp1);
			Z2.push_back(temp2);
			n += 1;
		}
		return;
	}

	void generateSs() {
		generateZs();
		double dt = op1.T / (double)M;
		int n = S1.size();
		while (n < N) {
			vector<double> temp1, temp2;
			double s1 = op1.S0, s2 = op2.S0;
			for (int i = 0;i < M;i++) {
				s1 *= exp((op1.r - op1.q - 0.5*op1.sigma*op1.sigma)*dt 
					+ op1.sigma*sqrt(dt)*Z1[n][i]);
				s2 *= exp((op2.r - op2.q - 0.5*op2.sigma*op2.sigma)*dt
					+ op2.sigma*sqrt(dt)*(rho*Z1[n][i]+sqrt(1-rho*rho)*Z2[n][i]));
				temp1.push_back(s1);
				temp2.push_back(s2);
			}
			S1.push_back(temp1);
			S2.push_back(temp2);
			n += 1;
		}
		return;
	}

	void generateVs() {
		generateSs();
		int n = Vs.size();
		while (n < N) {
			double max_value = op1.S0+op2.S0;
			for (int i = 0;i < M;i++) {
				max_value = max(max_value, S1[n][i] + S2[n][i]);
			}
			double v = exp(-op1.r*op1.T)*max(max_value - op1.K, 0.0);
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
