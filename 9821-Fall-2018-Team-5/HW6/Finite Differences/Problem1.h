#pragma once
#include"FiniteDifferencesPricer.h"

void problem1() {
	Option op(European, put, 42, 40, 0.5, 0.3, 0.02, 0.04);
	FiniteDifferencePircer pricer(op, 0.45, 4);
	
	ofstream out;
	out.open("1_data.csv");
	out << setprecision(10) << op.BSvalue << "," << op.BSdelta << endl << endl;

	for (int m = 0;m <= pricer.M;++m) {
		for (int n = 0;n <= pricer.N;++n) {
			out << setprecision(10) << pricer.Domin[n][m] << ",";
		}
		out << endl;
	}
	out << endl << endl;

	for (int M_ = 4;M_ <= 256;M_ *= 4) {
		FiniteDifferencePircer p(op, 0.45, M_);
		p.computeRMS();
		p.computeGreeks();
		out << setprecision(10) << p.error_pointwise << "," << "," <<
			p.error_pointwise2 << "," << "," <<
			p.RMS << "," << "," << 
			p.delta_approx << "," << p.gamma_approx << "," << p.theta_approx << endl;
	}
	out.close();
	
	return;
}
