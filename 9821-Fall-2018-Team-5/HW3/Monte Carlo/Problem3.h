#pragma once
#include"PathDependentcpp.h"
#include"NonPathDependent.h"
#include"NonePathDependentForProblem3.h"
#include"Generator.h"

void T3() {
	Option p(European, put, 60, 65, 0.5, 0.2, 0.0, 0.03);
	shared_ptr<Generator> gen_ptr1 = make_shared<InverseTransform>();
	shared_ptr<Generator> gen_ptr2 = make_shared<AcceptReject>();
	shared_ptr<Generator> gen_ptr3 = make_shared<BoxMuller>();
	int N = 10000;
	NonPathDependentForProblem3 pricer_1(N, gen_ptr1, p);
	NonPathDependentForProblem3 pricer_2(N, gen_ptr2, p);
	NonPathDependentForProblem3 pricer_3(N, gen_ptr3, p);
	
	std::ofstream out;
	out.open("3_data.csv");
	out << setprecision(10) << "BS value," << p.BSvalue << endl << endl;
	
	dynamic_pointer_cast<InverseTransform, Generator>(pricer_1.gen)->U.seed = 1;
	dynamic_pointer_cast<InverseTransform, Generator>(pricer_1.gen)->count = 0;
	for (N = 10000;N <= 5120000;N *= 2) {
		cout << N << endl;
		pricer_1.N = N;
		auto temp = pricer_1.MCvalue();
		out << N << "," << temp << "," << abs(p.BSvalue - temp) << endl;
	}
	out << endl;
	
	dynamic_pointer_cast<AcceptReject, Generator>(pricer_2.gen)->U.seed = 1;
	dynamic_pointer_cast<AcceptReject, Generator>(pricer_2.gen)->count = 0;
	for (N = 10000;N <= 5120000;N *= 2) {
		cout << N << endl;
		pricer_2.N = N;
		auto temp = pricer_2.MCvalue();
		out << N << ","<< pricer_2.Vs.size()<<"," << temp << "," << abs(p.BSvalue - temp) << endl;
	}
	out << endl;
	
	dynamic_pointer_cast<BoxMuller, Generator>(pricer_3.gen)->U.seed = 1;
	dynamic_pointer_cast<BoxMuller, Generator>(pricer_3.gen)->count = 0;
	for (N = 10000;N <= 5120000;N *= 2) {
		cout << N << endl;
		pricer_3.N = N;
		auto temp = pricer_3.MCvalue();
		out << N << "," << pricer_3.Vs.size() << "," << temp << "," << abs(p.BSvalue - temp) << endl;
	}

	out.close();

	
	return;
}