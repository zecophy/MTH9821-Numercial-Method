#pragma once
#include"Generator.h"

class NonPathDependentForProblem3 {
public:
	int N;
	shared_ptr<Generator> gen;
	vector<double> Zs;
	vector<double> Ss;
	vector<double> Vs;
	Option op;

	NonPathDependentForProblem3(int N_, shared_ptr<Generator> gen_, Option op_) :N(N_), gen(gen_), op(op_) {}

	void generateZs() {
		int n = Zs.size();
		while (1) {
			// generate one more
			if (gen->countNum() > N)break;
			Zs.push_back(gen->operator()());
			n += 1;
			

		} 
		return;
	}

	void generateSs() {
		generateZs();
		int n = Ss.size();
		// because generate one more
		while (n<Zs.size()-1) {
			double z = Zs[n];
			auto temp = op.S0*exp((op.r - op.q - op.sigma*op.sigma / 2)*op.T + op.sigma*sqrt(op.T)*z);
			Ss.push_back(temp);
			n += 1;
		}
		return;
	}

	void generateVs() {
		generateSs();
		int n = Vs.size();
		while (n<Ss.size()) {
			double s = Ss[n];
			Vs.push_back(exp(-op.r*op.T)*op.payoff(s));
			n += 1;
		}
		return;
	}

	double MCvalue() {
		generateVs();
		return accumulate(Vs.begin(), Vs.end(), 0.0) / (double)Vs.size();
	}

	vector<double> Deltas;
	vector<double> Vegas;

	void generateDeltas() {
		generateSs();
		int n = Deltas.size();
		while (n<N) {
			double s = Ss[n];
			Deltas.push_back(op.deltaMC(s));
			n += 1;
		}
		return;
	}

	void generateVegas() {
		generateSs();
		int n = Vegas.size();
		while (n<N) {
			double s = Ss[n];
			Vegas.push_back(op.vegaMC(s));
			n += 1;
		}
		return;
	}

	double MCdelta() {
		generateDeltas();
		return accumulate(Deltas.begin(), Deltas.end(), 0.0) / (double)N;
	}

	double MCvega() {
		generateVegas();
		return accumulate(Vegas.begin(), Vegas.end(), 0.0) / (double)N;
	}


};
