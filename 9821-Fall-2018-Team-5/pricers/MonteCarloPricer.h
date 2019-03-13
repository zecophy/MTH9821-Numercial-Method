/*
 * Copyright @ Chenyu Zhao
 * 
 * Supported　options: 1.European(Non-pathDepend) 2.Down-and-Out(PathDepend) 
 * 
 * Variance Reduction methods (Only for Non-pathDepend): 1.control variate 2.antithetic variable
 * 							3.moment matching 4.simultaneous moment matching and control variate
 * 
 * Use binary search to calculate implied volatility after have the MC value
 * 
 */ 

#pragma once
#include "Generator.h"

class MonteCarloPricer {
public:
	int N;
	int n,m;
	shared_ptr<Generator> gen;
	vector<double> Zs;
	vector<vector<double>> Ss;
	vector<double> Vs;
	Option op;
	
	double MCvalue;
	double MCdelta;
	double MCvega;
	
	double value_cv;// control variate
	double value_av;// antithetic variable
	double value_mm;// moment matching
	double value_sim;// simultaneous moment matching and control variate
	
	double imlVol;
	
	MonteCarloPricer(Option op_,int N_,shared_ptr<Generator> gen_):gen(gen_),op(op_){
		Reset(N_);
	}
	
	// reset N
	void Reset(int N_){
		N=N_;
		if(op.type==European){
			m=1;
			n=N;
		}
		if(op.type==DaO){// you can set how to choose m and n given N here
			m = 200;
			n = N / m;
		}
		
		PricerGo();
		return;
	}
	
	// start simulation
	void PricerGo(){
		generateVs();
		MCvalue=accumulate(Vs.begin(), Vs.end(), 0.0)/(double)n;
		
		ComputeGreeks();
		
		ControlVariate();
		AntitheticVariable();
		MomentMatching();
		SimultaneousMMandCV();
		
		//BSimpVol();
		return;
	}


/*
 * calculate implied vol
 */	
	void BSimpVol(double tol=0.000000001) {
		double value=MCvalue;
		auto op_=op;
		auto Cbs = [&](double sigma_)->double {
			op_.sigma = sigma_;
			op_.calculateBS();
			return op_.BSvalue;
		};

		double left = 0.0, right = 1;
		double fleft = Cbs(left) - value, fright = Cbs(right) - value;
		while (right - left >= tol) {
			double mid = (left + right) / 2;
			if ((Cbs(mid) - value) * fleft > 0) {
				left = mid;
			}
			else right = mid;
		}

		imlVol = (left + right) / 2;
		return;
	}

	
/*
 * Variance Reduction for NON-PATHDEPEND option
 */
	
	// for European options
	void ControlVariate() { 
		if(op.type!=European)return;
		
		// calculate b hat
		double S_mean = 0.0;
		for(auto path:Ss){
			S_mean+=path[m-1];
		}
		S_mean/=(double)n;
		double V_mean = accumulate(Vs.begin(), Vs.end(), 0.0) / (double)n;
		
		double nominator = 0.0;
		double dominator = 0.0;
		for (int i = 0;i < n;i++) {
			nominator += (Ss[i][m-1] - S_mean)*(Vs[i] - V_mean);
			dominator += (Ss[i][m-1] - S_mean)*(Ss[i][m-1] - S_mean);
		}
		double b_hat = nominator / dominator;

		// calculate Wi
		vector<double> Ws;
		int n_now = 0;

		double S_expectation = exp(op.r*op.T)*op.S0;
		while (n_now < n) {
			double w = Vs[n_now ] - b_hat * (Ss[n_now ][m-1] - S_expectation);
			Ws.push_back(w);
			n_now += 1;
		}

		value_cv = accumulate(Ws.begin(), Ws.end(), 0.0) / (double)n;

		return;
	}

	
	// for European Options
	void AntitheticVariable() {
		if(op.type!=European)return;
		
		// generate all MC values
		vector<double> Vs_anti;
		int n_now=0;
		while (n_now < n) {
			double s = op.S0*exp((op.r - op.q - 0.5*op.sigma*op.sigma)*op.T + op.sigma*sqrt(op.T)*(-Zs[n_now]));
			double v = op.payoff(s)*exp(-op.r*op.T);
			Vs_anti.push_back(v);
			n_now  += 1;
		}

		double sum = 0;
		for (int i = 0;i < n;i++) {
			sum += (Vs[i] + Vs_anti[i]) / 2.0;
		}
		value_av = sum / (double)n;

		return;
	}
	
	vector<double> Vs_mm;
	vector<double> Smm;
	// for European Options
	void MomentMatching() {
		if(op.type!=European)return;
		// calculate V^tilde
		double S_mean = 0.0;
		for(auto path:Ss){
			S_mean+=path[m-1];
		}
		S_mean/=(double)n;
		double S_expectation = exp(op.r*op.T)*op.S0;

		int n_now = 0;
		Vs_mm.clear();
		Smm.clear();

		while (n_now < n) {
			double s = Ss[n_now][m-1] * S_expectation / S_mean;
			Smm.push_back(s);
			double v = op.payoff(s)*exp(-op.r*op.T);
			Vs_mm.push_back(v);
			n_now += 1;
		}

		// calculate Vmm
		value_mm = accumulate(Vs_mm.begin(), Vs_mm.end(), 0.0) / (double)n;

		return;
	}

	// for European Options
	void SimultaneousMMandCV() {
		if(op.type!=European)return;
		double V_mean = accumulate(Vs_mm.begin(), Vs_mm.end(), 0.0) / (double)n;
		double S_expectation = exp(op.r*op.T)*op.S0;

		// calculate b hat
		double nominator = 0;
		double dominator = 0;
		for (int i = 0;i < n;i++) {
			nominator += (Smm[i] - S_expectation)*(Vs_mm[i] - V_mean);
			dominator += (Smm[i] - S_expectation)*(Smm[i] - S_expectation);
		}
		double b_hat = nominator / dominator;

		// calculate Wi
		vector<double> Ws_sim;
		int n_now = 0;
		while (n_now < n) {
			double w = Vs_mm[n_now] - b_hat * (Smm[n_now] - S_expectation);
			Ws_sim.push_back(w);
			n_now += 1;
		}

		// calculate Vsim
		value_sim = accumulate(Ws_sim.begin(), Ws_sim.end(), 0.0) / (double)n;

		return;
	}



private:
/*
 * Generating Function: you can modify when you have new options
 */
	void generateZs() {
		int t = Zs.size();
		while(t<N){
			Zs.push_back(gen->operator()());
			t += 1;
		}

		return;
	}

	void generateSs() {
		generateZs();
		
		Ss.clear();
		double dt = op.T / (double)m;
		
		int n_now=0;
		while (n_now < n) {
			vector<double> path;
			double s = op.S0;
			for (int j = 0;j < m;j++) {
				double z = Zs[n_now*m + j];
				s *= exp((op.r - op.q - op.sigma*op.sigma / 2)*dt + op.sigma*sqrt(dt)*z);
				path.push_back(s);
			}
			Ss.push_back(path);
			n_now++;
		}
	
		return;
	}

	void generateVs() {
		generateSs();
		Vs.clear();
		
		for (int i = 0;i < n;i++) {
			if(op.type==European){
				Vs.push_back(exp(-op.r*op.T)*op.payoff(Ss[i][m - 1]));
			}
			if(op.type==DaO){
				double min_S = *min_element(Ss[i].begin(), Ss[i].end());
				if (min_S < op.B)
					Vs.push_back(0.0);
				else
					Vs.push_back(exp(-op.r*op.T)*op.payoff(Ss[i][m - 1]));
			}
			
		}
		return;
	}


/*
 * Generating Greeks
 */
 
	void ComputeGreeks(){
		// delta
		generateDeltas();
		MCdelta=accumulate(Deltas.begin(), Deltas.end(), 0.0) / (double)n;
		
		// vega
		generateVegas();
		MCvega= accumulate(Vegas.begin(), Vegas.end(), 0.0) / (double)n;
		
		return;
	}
	
	
	vector<double> Deltas;
	vector<double> Vegas;

	void generateDeltas() {
		generateSs();
		int t = Deltas.size();
		while (t<n) {
			double s = Ss[t][m-1];
			Deltas.push_back(deltaMC(s));
			t += 1;
		}
		return;
	}

	void generateVegas() {
		generateSs();
		int t = Vegas.size();
		while (t<n) {
			double s = Ss[t][m-1];
			Vegas.push_back(vegaMC(s));
			t += 1;
		}
		return;
	}
	
	// used to calculate greeks given stock price for each simulation
	double deltaMC(double S) {
		if (op.cp == put) {
			if (op.K > S)
				return -exp(-op.r * op.T)*S / op.S0;
			else
				return 0.0;
		}
		if (op.cp == call) {
			if (S > op.K)
				return exp(-op.r * op.T)*S / op.S0;
			else
				return 0.0;
		}
		return -1;
	}

	double vegaMC(double S) {
		double z = (log(S / op.S0) - (op.r - op.q - op.sigma * op.sigma / 2)*op.T)
				/ op.sigma / sqrt(op.T);
		if (op.cp == put) {
			if (op.K > S)
				return -S * exp(-op.r * op.T)*(-op.sigma * op.T + sqrt(op.T)*z);
			else
				return 0;
		}
		if (op.cp == call) {
			if (S > op.K)
				return S * exp(-op.r * op.T)*(-op.sigma * op.T + sqrt(op.T)*z);
			else
				return 0;
		}
		return -1;
	}
};
