#pragma once
#include"NonPathDependent.h"

class PathDependent {
public:
	int N;
	int n;
	int m;
	shared_ptr<Generator> gen;
	vector<double> Zs;
	vector<vector<double>> Ss;
	vector<double> Vs;
	Option op;

	PathDependent(int N_, shared_ptr<Generator> gen_, Option op_) 
		:N(N_), gen(gen_), op(op_) {
		m = 200;n = N / m;
	}

	void generateZs() {
		int t = Zs.size();
		while (t<N) {
			Zs.push_back(gen->operator()());
			t += 1;
		}
		return;
	}

	void generateSs() {
		generateZs();
		int n_now = Ss.size();
		if (n_now != 0) {
			if (Ss[0].size() != m) {
				Ss.clear();
				n_now = 0;
			}
		}
		
		
		double dt = op.T / (double)m;

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
		//only for DoA
		{	
			Vs.clear();
			for (int i = 0;i < n;i++) {
				double min_S = *min_element(Ss[i].begin(), Ss[i].end());
				if (min_S < op.B)Vs.push_back(0.0);
				else
					Vs.push_back(exp(-op.r*op.T)*op.payoff(Ss[i][m - 1]));
			}
		}
		return;
	}

	double MCvalue() {
		generateVs();
		return accumulate(Vs.begin(), Vs.end(), 0.0) / (double)Vs.size();
	}

	

};
