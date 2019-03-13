#include"FiniteDifferencesPricer.h"

int main() {
	double exact_value = 3.880035021683392;
	Option op(American, put, 42, 40, 1, 0.32, 0.02, 0.04);
	FiniteDifferencePircer pricer(op, 0.45, 4, CN_SOR, exact_value);

	ofstream out;
	out.open("1_data.csv");
	out << setprecision(10) << pricer.v_approx << "," << pricer.v_approx_red << endl << endl;

	for (int m = 0;m <= pricer.M;++m) {
		for (int n = 0;n <= pricer.N;++n) {
			out << setprecision(10) << pricer.Domin[n][m] << ",";
		}
		out << endl;
	}
	out << endl << endl;

	for (int M_ = 4;M_ <= 256;M_ *= 4) {
		cout << M_ << endl;
		FiniteDifferencePircer p(op, 0.45, M_, CN_SOR, exact_value);
		p.computeGreeks();
		p.computeRed();
		out << setprecision(10) << p.error_pointwise << "," << "," <<
			p.error_pointwise2 << "," << "," <<
			p.delta_approx << "," << p.gamma_approx << "," << p.theta_approx << "," <<
			p.v_approx_red << "," << p.error_pointwise_red << endl;
	}
	out << endl << endl;

	for (int M_ = 4;M_ <= 256;M_ *= 4) {
		cout << M_ << endl;
		FiniteDifferencePircer p(op, 5.0, M_, CN_SOR, exact_value);
		p.computeGreeks();
		p.computeRed();
		out << setprecision(10) << p.error_pointwise << "," << "," <<
			p.error_pointwise2 << "," << "," <<
			p.delta_approx << "," << p.gamma_approx << "," << p.theta_approx << "," <<
			p.v_approx_red << "," << p.error_pointwise_red << endl;
	}
	out << endl << endl;

	out.close();
	system("pause");
	return 0;
}