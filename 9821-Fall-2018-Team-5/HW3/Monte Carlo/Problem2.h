#pragma once
#include"PathDependentcpp.h"

void T2() {
	Option dao(DaO, call, 40, 40, 0.75, 0.3, 0.01, 0.02, 35);
	//dao.calculateBS();
	int N = 10000;
	shared_ptr<Generator> gen_ptr1 = make_shared<InverseTransform>();
	shared_ptr<Generator> gen_ptr2 = make_shared<InverseTransform>();
	PathDependent pricer_1(N, gen_ptr1, dao);
	PathDependent pricer_2(N, gen_ptr2, dao);

	std::ofstream out;
	out.open("2_data.csv");
	out << setprecision(10) << "BS value," << dao.BSvalue << endl << endl;

	/*
	pricer_1.generateVs();
	for (auto item : pricer_1.Vs)
		cout << item << endl;
		*/
	dynamic_pointer_cast<InverseTransform, Generator>(pricer_1.gen)->U.seed = 1;
	for (;N <= 5120000;N *= 2) {
		cout << N << endl;
		pricer_1.N = N;
		pricer_1.m = 200;
		pricer_1.n = N / 200;
		double temp = pricer_1.MCvalue();
		out << setprecision(10) << pricer_1.m << "," << pricer_1.n << ","
			<< temp <<","<< abs(pricer_1.op.BSvalue - temp) << endl;
	}


	dynamic_pointer_cast<InverseTransform, Generator>(pricer_2.gen)->U.seed = 1;
	out << endl<<endl;
	for (N=10000;N <= 5120000;N *= 2) {
		cout << N << endl;
		int mk = ceil(pow((double)N, 1.0 / 3.0)*pow((double)0.75, 2.0 / 3.0));
		int nk = floor((double)N / (double)mk);
		pricer_2.N = N;
		pricer_2.m = mk;
		pricer_2.n = nk;
		double temp = pricer_2.MCvalue();
		out << setprecision(10) << pricer_2.m << "," << pricer_2.n << ","
			<< temp << "," << abs(pricer_2.op.BSvalue - temp) << endl;
	}
	

	out.close();


	//system("pause");
	return;
}
