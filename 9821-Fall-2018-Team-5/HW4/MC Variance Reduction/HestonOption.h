#pragma once
#include"Generator.h"

class HestonOption:public Option {
public:
	double v0;
	double lambda;
	double eta;
	double v_mean;
	double rho;

	HestonOption(Type type_,CP cp_,double S0_,double K_,double T_,double q_,double r_,
		double v0_,double lambda_,double v_mean_,double eta_,double rho_): Option(type_,cp_,
		S0_,K_,T_,0,q_,r_),lambda(lambda_),eta(eta_),v_mean(v_mean_),rho(rho_),v0(v0_){}
};

class HestonOptionPricer {
public:
	BoxMuller gen;
	int N, m;
	HestonOption op;

	vector<vector<double>> Z1;
	vector<vector<double>> Z2;
	vector<vector<double>> S;
	vector<vector<double>> Vol;
	vector<double> Vs;

	double value;
	double imlVol;

	HestonOptionPricer(HestonOption op_,int m_):op(op_),m(m_){}

	void generateZs() {
		int n = Z1.size();
		while (n < N) {
			vector<double> temp1;
			vector<double> temp2;
			for (int i = 0;i < m;i++) {
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
		double dt = op.T / (double)m;
		int n = S.size();
		while (n < N) {
			vector<double> temp1, temp2;
			double s = op.S0, v = op.v0;
			for (int i = 0;i < m;i++) {
				double vplus = max(v, 0.0);
				v = vplus - op.lambda*(vplus - op.v_mean)*dt + op.eta*sqrt(vplus*dt)*
					(op.rho*Z1[n][i] + sqrt(1 - op.rho*op.rho)*Z2[n][i]);
				s *= exp((op.r - op.q - 0.5*vplus)*dt+ sqrt(vplus*dt)*Z1[n][i] );
				temp1.push_back(s);
				temp2.push_back(v);
			}
			S.push_back(temp1);
			Vol.push_back(temp2);
			n += 1;
		}
		return;
	}

	void generateVs() {
		generateSs();
		int n = Vs.size();
		while (n < N) {
			double v = exp(-op.T*op.r)*op.payoff(S[n][m - 1]);
			Vs.push_back(v);
			n += 1;
		}
		return;
	}

	void BSimpVol() {
		auto Cbs = [=](double sigma_)->double {
			op.sigma = sigma_;
			op.calculateBS();
			return op.BSvalue;
		};

		double left = 0.0, right = 1;
		double fleft = Cbs(left) - value, fright = Cbs(right) - value;
		while (right - left >= 0.000000001) {
			double mid = (left + right) / 2;
			if ((Cbs(mid) - value) * fleft > 0) {
				left = mid;
			}
			else right = mid;
		}

		imlVol = (left + right) / 2;
		return;
	}

	void PricerGo(int N_) {
		N = N_;
		generateVs();
		value = accumulate(Vs.begin(), Vs.end(), 0.0) / (double)N;
		BSimpVol();
		return;
	}

};
