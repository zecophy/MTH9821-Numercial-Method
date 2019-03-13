#pragma once
#include"FiniteDifferencesPricer.h"

void problem2() {
	double exact_value = 2.380407113545689;
	Option op(American, put, 42, 40, 0.5, 0.3, 0.02, 0.04);
	FiniteDifferencePircer pricer(op, 0.45, 4, exact_value);

	ofstream out;
	out.open("2_data.csv");
	out << setprecision(10) << pricer.v_approx << "," << pricer.v_approx_red << endl << endl;

	for (int m = 0;m <= pricer.M;++m) {
		for (int n = 0;n <= pricer.N;++n) {
			out << setprecision(10) << pricer.Domin[n][m] << ",";
		}
		out << endl;
	}
	out << endl << endl;

	/*
	for (int m = 0;m <= pricer.M;++m) {
		for (int n = 0;n <= pricer.N;++n) {
			out << setprecision(10) << pricer.early_ex[n][m] << ",";
		}
		out << endl;
	}
	out << endl << endl;
	*/

	for (int M_ = 4;M_ <= 256;M_ *= 4) {
		FiniteDifferencePircer p(op, 0.45, M_, exact_value);
		//p.computeRMS();
		p.computeGreeks();
		p.computeRed();
		out << setprecision(10) << p.error_pointwise << "," << "," <<
			p.error_pointwise2 << "," << "," <<
			p.delta_approx << "," << p.gamma_approx << "," << p.theta_approx << "," <<
			p.v_approx_red << "," << p.error_pointwise_red << endl;
	}
	out << endl << endl;


	FiniteDifferencePircer pp(op, 0.45, 16, exact_value);
	pp.computeEarly();
	for (int i = 0;i <= pp.M;++i) {
		out << setprecision(10) << i << "," << pp.early_t[i] << "," << pp.early_s[i] << endl;
	}


	out.close();

	return;
}
