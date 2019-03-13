#pragma once
#include"FiniteDifferencesPricer.h"

void problem1() {
	Option op(European, put, 42, 40, 0.5, 0.3, 0.02, 0.04);
	FiniteDifferencePircer pricer_1(op, 0.45, 4, Backward);
	FiniteDifferencePircer pricer_2(op, 0.45, 4, CN);
	
	ofstream out;
	out.open("1_data.csv");

	for (auto & pricer : { pricer_1,pricer_2 }) {
		if (pricer.fd_type == 1)out << "Backward" << endl;
		else out << "CN" << endl;
		for (int m = 0;m <= pricer.M;++m) {
			for (int n = 0;n <= pricer.N;++n) {
				out << setprecision(10) << pricer.Domin[n][m] << ",";
			}
			out << endl;
		}
		out << endl << endl;
	}

	// the first part: output domin end
	///////////////////////////////////////////
	

	for (double alpha : {0.45, 5.0}) {
		for (auto fd_t : { Backward,CN }) {
			out << alpha << ",";
			if (fd_t == 1)out << "Backward" << endl;
			else out << "CN" << endl;
			for (int M_ = 4;M_ <= 256;M_ *= 4) {
				FiniteDifferencePircer p(op, alpha, M_, fd_t);
				p.computeRMS();
				p.computeGreeks();
				out << setprecision(10) << p.error_pointwise << "," << "," <<
					p.error_pointwise2 << "," << "," <<
					p.RMS << "," << "," <<
					p.delta_approx << "," << p.gamma_approx << "," << p.theta_approx << endl;
			}
			out << endl << endl << endl;
		}
	}
	
	
	out.close();
	
	return;
}
