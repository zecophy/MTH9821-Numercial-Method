#pragma once
#include"EuropeanOption.h"

void problem1() {
	Option op(European, put, 50, 55, 0.75, 0.25, 0.0, 0.03);
	EuropeanOptionPricer pricer(op);

	ofstream out;
	out.open("1_data.csv");
	out << op.BSvalue << endl << endl;

	//5120000
	for (int N = 10000;N <= 5120000;N *= 2) {
		cout << N << endl;
		pricer.VarianceReductionGo(N);
		double Vbs = op.BSvalue;
		out << setprecision(10) << N << ","
			<< pricer.Vav << "," << abs(pricer.Vav - Vbs) << ","
			<< pricer.Vsim << "," << abs(pricer.Vsim - Vbs) << ","
			<< pricer.Vcv << "," << abs(pricer.Vcv - Vbs) << ","
			<< pricer.Vmm << "," << abs(pricer.Vmm - Vbs) << endl;
	}

	out.close();

	return;
}
