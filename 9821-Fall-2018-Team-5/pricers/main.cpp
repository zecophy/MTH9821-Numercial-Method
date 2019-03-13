#include <iostream>
#include <vector>
#include <memory>
#include "FiniteDifferencesPricer.h"
#include "DividendFiniteDifferencePricer.h"
#include "MonteCarloPricer.h"
#include "TreePricer.h"
using namespace std;

int main(int argc, char **argv) {
	Option op1(European,call,40,40,1,0.25,0.01,0.03);
	Option op2(American,put,40,40,1,0.25,0.01,0.03);
	Option op3(DaO,call,40,40,1,0.25,0.01,0.03, 30);
	
	
	/*
	int N=50;
	FiniteDifferencePircer pricer11(op2,0.45,N,Forward);
	FiniteDifferencePircer pricer12(op2,0.45,N,Backward);
	FiniteDifferencePircer pricer13(op2,0.45,N,CN);
	FiniteDifferencePircer pricer14(op2,0.45,N,CN_SOR);
	
	vector<FiniteDifferencePircer> pricers1{pricer11,pricer12,pricer13,pricer14};
	
	for(auto& pricer:pricers1){
		for(int m=0;m<=pricer.M;++m){
			for(int n=0;n<=pricer.N;++n){
				cout<<pricer.Domin[n][m]<<", ";
			}
			cout<<endl;
		}
		cout<<"\n*********************************\n";
	}
	
	cout<<"\n\n\n";
	cout<<pricers1[0].op.BSvalue<<endl;
	for(auto& pricer:pricers1){
		cout<<pricer.v_approx<<", "<<pricer.error_pointwise<<", "<<pricer.v_approx2<<", "<<pricer.error_pointwise2<<endl;
	}
	cout<<"\n\n\n";
	for(auto& pricer:pricers1){
		cout<<pricer.RMS<<endl;
	}
	
	cout<<"\n\n\n";
	cout<<pricers1[0].op.BSdelta<<", "<<pricers1[0].op.BSgamma<<", "<<pricers1[0].op.BStheta<<endl;
	for(auto& pricer:pricers1){
		cout<<pricer.delta_approx<<", "<<pricer.gamma_approx<<", "<<pricer.theta_approx<<endl;
	}
	
	
	cout<<"\n\n\n";
	for(auto& pricer:pricers1){
		cout<<pricer.v_approx_red<<endl;
	}
	
	cout<<"\n\n\n";
	for(auto& pricer:pricers1){
		for(int i=0;i<pricer.early_t.size();++i){
			cout<<pricer.early_t[i]<<", "<<pricer.early_s[i]<<endl;
		}
		cout<<"\n\n";
	}
	*/
	
	
	
	
	
	/*
	shared_ptr<Generator> gen1=make_shared<BoxMuller>();
	shared_ptr<Generator> gen2=make_shared<BoxMuller>();
	MonteCarloPricer pricer1(op1,500000,gen1);
	MonteCarloPricer pricer2(op3,50000,gen2);
	
	vector<MonteCarloPricer> p{pricer1};
	

	for(auto& item : p){
		item.PricerGo();cout<<item.op.BSvalue<<endl;
		cout<<item.MCvalue<<", "<<item.MCdelta<<", "<<item.MCvega<<endl;
		cout<<item.value_cv <<", "<<item.value_av<<", "<< item.value_mm<<", "<< item.value_sim<<endl;
	}
	*/
	/*
	Option dao(DaO, call, 40, 40, 0.75, 0.3, 0.01, 0.02, 35);
	int N = 10000;
	shared_ptr<Generator> gen_ptr1 = make_shared<InverseTransform>();
	shared_ptr<Generator> gen_ptr2 = make_shared<InverseTransform>();
	MonteCarloPricer pricer_1(dao,N, gen_ptr1);
	MonteCarloPricer pricer_2(dao,N, gen_ptr2);

	std::ofstream out;
	out.open("2_data.csv");
	out << setprecision(10) << "BS value," << dao.BSvalue << endl << endl;

	/*
	pricer_1.generateVs();
	for (auto item : pricer_1.Vs)
		cout << item << endl;
		*/
		/*
	dynamic_pointer_cast<InverseTransform, Generator>(pricer_1.gen)->U.seed = 1;
	for (;N <= 5120000;N *= 2) {
		cout << N << endl;
		pricer_1.Reset(N);
		pricer_1.PricerGo();

		double temp = pricer_1.MCvalue;
		out << setprecision(10) << pricer_1.m << "," << pricer_1.n << ","
			<< temp <<","<< abs(pricer_1.op.BSvalue - temp) << endl;
	}




	out.close();
	*/
	/*
	auto op=op1;
	BinomialTreePricer pricer(op,1000);
	cout<<op.BSvalue<<", "<<op.BSdelta<<", "<<op.BSgamma<<", "<<op.BStheta<<endl;
	pricer.Binomial();
	cout<<pricer.Treevalue<<", "<<pricer.Treedelta<<", "<<pricer.Treegamma<<", "<<pricer.Treetheta<<endl;
	pricer.AverageBinomial();
	cout<<pricer.Treevalue<<", "<<pricer.Treedelta<<", "<<pricer.Treegamma<<", "<<pricer.Treetheta<<endl;
	pricer.BSBinomial();
	cout<<pricer.Treevalue<<", "<<pricer.Treedelta<<", "<<pricer.Treegamma<<", "<<pricer.Treetheta<<endl;
	pricer.BBSR();
	cout<<pricer.Treevalue<<", "<<pricer.Treedelta<<", "<<pricer.Treegamma<<", "<<pricer.Treetheta<<endl;
	*/
	auto op=op1;
	TrinomialTreePricer pricer(op,1000);
	cout<<op.BSvalue<<", "<<op.BSdelta<<", "<<op.BSgamma<<", "<<op.BStheta<<endl;
	pricer.Trinomial();
	cout<<pricer.Treevalue<<", "<<pricer.Treedelta<<", "<<pricer.Treegamma<<", "<<pricer.Treetheta<<endl;
	pricer.BSTrinomial();
	cout<<pricer.Treevalue<<", "<<pricer.Treedelta<<", "<<pricer.Treegamma<<", "<<pricer.Treetheta<<endl;
	pricer.TBSR();
	cout<<pricer.Treevalue<<", "<<pricer.Treedelta<<", "<<pricer.Treegamma<<", "<<pricer.Treetheta<<endl;
	return 0;
}
