#pragma once
#include"Generator.h"

class EuropeanOptionPricer {
public:
	int N;
	Option op;
	BoxMuller gen;
	vector<double> Zs;
	vector<double> Zs_minus;
	vector<double> Ss;
	vector<double> Ss_minus;
	vector<double> Smm;
	vector<double> Vs;
	vector<double> Vs_minus;
	vector<double> Vs_mm;
	vector<double> Ws;
	vector<double> Ws_sim;

	double Vcv;
	double Vav;
	double Vmm;
	double Vsim;
	
	EuropeanOptionPricer(Option op_):op(op_){}

	void refresh() {
		Zs.clear();
		Zs_minus.clear();
		Ss.clear();
		Ss_minus.clear();
		Smm.clear();
		Vs.clear();
		Vs_minus.clear();
		Ws.clear();
		Ws_sim.clear();
		gen.refesh();
	}

	void generateZs() {
		int n = Zs.size();
		while (n < N) {
			double z = gen();
			Zs.push_back(z);
			Zs_minus.push_back(-z);
			n += 1;
		}
		return;
	}

	//must be used firstly
	void ControlVariate() { 
		// generate all MC values
		generateZs();
		int n = Ss.size();
		while (n < N) {
			double s = op.S0*exp((op.r - op.q - 0.5*op.sigma*op.sigma)*op.T + op.sigma*sqrt(op.T)*Zs[n]);
			Ss.push_back(s);
			double v = op.payoff(s)*exp(-op.r*op.T);
			Vs.push_back(v);
			n += 1;
		}

		// calculate b hat
		double S_mean = accumulate(Ss.begin(), Ss.end(), 0.0) / (double)N;
		double V_mean = accumulate(Vs.begin(), Vs.end(), 0.0) / (double)N;
		
		double nominator = 0;
		double dominator = 0;
		for (int i = 0;i < N;i++) {
			nominator += (Ss[i] - S_mean)*(Vs[i] - V_mean);
			dominator += (Ss[i] - S_mean)*(Ss[i] - S_mean);
		}
		double b_hat = nominator / dominator;

		// calculate Wi
		// because we update S_mean with every new N, we must recalculate all S
		n = 0;
		Ws.clear();
		double S_expectation = exp(op.r*op.T)*op.S0;
		while (n < N) {
			double w = Vs[n] - b_hat * (Ss[n] - S_expectation);
			Ws.push_back(w);
			n += 1;
		}

		// calculate Vcv
		Vcv = accumulate(Ws.begin(), Ws.end(), 0.0) / (double)N;

		return;
	}

	//must be used secondly
	void AntitheticVariable() {
		// generate all MC values
		int n = Ss_minus.size();
		while (n < N) {
			double s = op.S0*exp((op.r - op.q - 0.5*op.sigma*op.sigma)*op.T + op.sigma*sqrt(op.T)*Zs_minus[n]);
			Ss_minus.push_back(s);
			double v = op.payoff(s)*exp(-op.r*op.T);
			Vs_minus.push_back(v);
			n += 1;
		}

		double sum = 0;
		for (int i = 0;i < N;i++) {
			sum += (Vs[i] + Vs_minus[i]) / 2.0;
		}
		Vav = sum / (double)N;

		return;
	}

	//must be used secondly
	void MomentMatching() {
		// calculate V^tilde
		double S_mean = accumulate(Ss.begin(), Ss.end(), 0.0) / (double)N;
		double S_expectation = exp(op.r*op.T)*op.S0;

		// because we update S_mean with every new N, we must recalculate all S
		int n = 0;
		Vs_mm.clear();
		Smm.clear();
		while (n < N) {
			double s = Ss[n] * S_expectation / S_mean;
			Smm.push_back(s);
			double v = op.payoff(s)*exp(-op.r*op.T);
			Vs_mm.push_back(v);
			n += 1;
		}

		// calculate Vmm
		Vmm = accumulate(Vs_mm.begin(), Vs_mm.end(), 0.0) / (double)N;

		return;
	}

	// must be used in the end
	void SimultaneousMMandCV() {
		double V_mean = accumulate(Vs_mm.begin(), Vs_mm.end(), 0.0) / (double)N;
		double S_expectation = exp(op.r*op.T)*op.S0;

		// calculate b hat
		double nominator = 0;
		double dominator = 0;
		for (int i = 0;i < N;i++) {
			nominator += (Smm[i] - S_expectation)*(Vs_mm[i] - V_mean);
			dominator += (Smm[i] - S_expectation)*(Smm[i] - S_expectation);
		}
		double b_hat = nominator / dominator;

		// calculate Wi
		// because we update S_mean with every new N, we must recalculate all S
		int n = 0;
		Ws_sim.clear();
		while (n < N) {
			double w = Vs_mm[n] - b_hat * (Smm[n] - S_expectation);
			Ws_sim.push_back(w);
			n += 1;
		}

		// calculate Vsim
		Vsim = accumulate(Ws_sim.begin(), Ws_sim.end(), 0.0) / (double)N;

		return;
	}

	void VarianceReductionGo(int N_) {
		N = N_;
		ControlVariate();
		AntitheticVariable();
		MomentMatching();
		SimultaneousMMandCV();
		return;
	}

	
};
