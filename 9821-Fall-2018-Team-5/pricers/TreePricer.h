/*
 * Copyright @ Chenyu Zhao
 * 
 * Supported options: 1.European) 2.American
 * 
 * Supported methods: 1. BinomialTree, AverageBinomialTree, BSBinomialTree, BBSR
 * 					  2. TrinomialTree, BSTrinomialTree, TBSR
 */ 

#pragma once
#include "Option.h"

class BinomialTreePricer{
public:
	Option op;
	int N;
	
	double dt;
	double disc;
	double u;
	double d;
	double pu;
	double pd;
	
	std::vector<std::vector<double>> store{ {0},{0,0},{0,0,0} };
	
	double Treevalue;
	double Treedelta,Treegamma,Treetheta;

	BinomialTreePricer(Option op_, int N_):op(op_),N(N_){
		dt = op.T / N;
		disc = std::exp(-op.r*dt);
		u = std::exp(op.sigma*std::sqrt(dt));
		d = 1 / u;
		pu = (std::exp((op.r - op.q)*dt) - d) / (u - d);
		pd = 1 - pu;
	}
	
	
	
	void Binomial(){
		// initial last period price
		double Smax = op.S0 * std::pow(u, N);
		std::vector<double> V(N + 1);
		for (int i = 0;i <= N;i++) {
			V[i]=op.payoff(Smax);
			Smax /= (u*u);
		}
		
		// if european option
		if (op.type == European){
			for (int k = N - 1;k >= 0;k--) {
				for (int i = 0;i <= k;i++) {
					V[i] = disc * (pu*V[i] + pd * V[i + 1]);
					if (k <= 2)store[k][i] = V[i];
				}
			}
		}

		// if american option
		else if (op.type == American) {
			for (int k = N - 1;k >= 0;k--) {
				double Smax2 = op.S0*std::pow(u, k);
				for (int i = 0;i <= k;i++) {
					V[i] = std::max(disc*(pu*V[i] + pd * V[i + 1]), op.payoff(Smax2));				

					if (k <= 2)store[k][i] = V[i];
					Smax2 /= (u*u);
				}
			}
		}
		
		Treevalue=V[0];
		Treedelta = (store[1][0] - store[1][1]) / (op.S0*u - op.S0*d);
		Treegamma = ((store[2][0] - store[2][1]) / (op.S0*u*u - op.S0)
					- (store[2][1] - store[2][2]) / (op.S0 - op.S0*d*d)) / ((op.S0*u*u - op.S0*d*d) / 2);
		Treetheta = (store[2][1] - store[0][0]) / 2 / dt;

		return;
	}
	
	
	void AverageBinomial(){
		BinomialTreePricer a(op,N);
		BinomialTreePricer b(op,N+1);
		a.Binomial();
		b.Binomial();
		Treevalue=(a.Treevalue+b.Treevalue)/2.0;
		Treedelta=(a.Treedelta+b.Treedelta)/2.0;
		Treegamma=(a.Treegamma+b.Treegamma)/2.0;
		Treetheta=(a.Treetheta+b.Treetheta)/2.0;
		return;
	}
	
	
	void BSBinomial(){
		double Smax = op.S0 * std::pow(u, N-1);
		std::vector<double> V(N);

		for (int i = 0;i < N;i++) {
			Option temp = op;
			temp.T = dt;
			temp.S0 = Smax;
			temp.calculateBS();
			if (op.type == European)
				V[i] = temp.BSvalue;
			if(op.type==American)
				V[i]=max(temp.BSvalue,op.payoff(Smax));
			Smax /= (u*u);
		}


		// if european option
		if (op.type == European) {
			for (int k = N - 2;k >= 0;k--) {
				for (int i = 0;i <= k;i++) {
					V[i] = disc * (pu*V[i] + pd * V[i + 1]);
					if (k <= 2)store[k][i] = V[i];
				}
			}
		}

		// if american option
		else if (op.type == American) {
			for (int k = N - 2;k >= 0;k--) {
				double Smax2 = op.S0*std::pow(u, k);
				for (int i = 0;i <= k;i++) {
					V[i] = std::max(disc*(pu*V[i] + pd * V[i + 1]), op.payoff(Smax2));
				
					if (k <= 2)store[k][i] = V[i];
					Smax2 /= (u*u);
				}
			}
		}
		
		Treevalue=V[0];
		Treedelta = (store[1][0] - store[1][1]) / (op.S0*u - op.S0*d);
		Treegamma = ((store[2][0] - store[2][1]) / (op.S0*u*u - op.S0)
					- (store[2][1] - store[2][2]) / (op.S0 - op.S0*d*d)) / ((op.S0*u*u - op.S0*d*d) / 2);
		Treetheta = (store[2][1] - store[0][0]) / 2 / dt;
		return;
		
	}
	
	void BBSR(){
		BinomialTreePricer a(op,N);
		BinomialTreePricer b(op,N/2);
		a.BSBinomial();
		b.BSBinomial();
		Treevalue=(2*a.Treevalue-b.Treevalue);
		Treedelta=(2*a.Treedelta-b.Treedelta);
		Treegamma=(2*a.Treegamma-b.Treegamma);
		Treetheta=(2*a.Treetheta-b.Treetheta);
		return;
	}
	
};





class TrinomialTreePricer{
public:
	Option op;
	int N;
	
	double dt;
	double disc;
	double u;
	double m;
	double d;
	double pu;
	double pm;
	double pd;
	
	std::vector<std::vector<double>> store{ { 0 },{ 0,0,0 },{ 0,0,0,0,0 } };
	
	double Treevalue;
	double Treedelta,Treegamma,Treetheta;
	
	
	TrinomialTreePricer(Option op_, int N_):op(op_),N(N_){
		dt = op.T / N;
		disc = std::exp(-op.r*dt);
		u = std::exp(op.sigma*std::sqrt(3*dt));
		m = 1;
		d = 1 / u;
		pu = 1.0 / 6.0 + (op.r - op.q - 0.5*op.sigma*op.sigma)*std::sqrt(dt / 12 / op.sigma / op.sigma);
		pm = 2.0 / 3.0;
		pd = 1 - pu - pm;
	}
	
	void Trinomial(){
		double Smax = op.S0 * std::pow(u, N);
		std::vector<double> V(2 * N + 1);
		for (int i = 0;i <= 2*N;i++) {
			V[i] = op.payoff(Smax);
			Smax /= u;
		}
		
		// if european option
		if (op.type == European) {
			for (int k = N - 1;k >= 0;k--) {
				for (int i = 0;i <= 2*k;i++) {
					V[i] = disc * (pu*V[i] + pm * V[i + 1] + pd * V[i + 2]);
					if (k <= 2)store[k][i] = V[i];
				}
			}
		}

		// if american option
		else if (op.type == American) {
			for (int k = N - 1;k >= 0;k--) {
				double Smax2 = op.S0*std::pow(u, k);
				for (int i = 0;i <= 2*k;i++) {
					V[i] = std::max(disc * (pu*V[i] + pm * V[i + 1] + pd * V[i + 2]), op.payoff(Smax2));

					if (k <= 2)store[k][i] = V[i];
					Smax2 /= u;
				}
			}
		}
		
		Treevalue=V[0];
		Treedelta = (store[1][0] - store[1][2]) / (op.S0*u - op.S0*d);
		Treegamma = ((store[2][0] - store[2][2]) / (op.S0*u*u - op.S0)
					- (store[2][2] - store[2][4]) / (op.S0 - op.S0*d*d)) / (op.S0*u - op.S0*d);
		Treetheta = (store[1][1] - store[0][0]) /  dt;
		return;
	}
	
	void BSTrinomial(){
		double Smax = op.S0 * std::pow(u, N - 1);
		std::vector<double> V(2*N-1);
		
		for (int i = 0;i < 2*N-1;i++) {
			Option temp = op;
			temp.T = dt;
			temp.S0 = Smax;
			temp.calculateBS();
			if (op.type == European)
				V[i] = temp.BSvalue;
			if(op.type==American)
				V[i]=max(temp.BSvalue,op.payoff(Smax));
			Smax /= u;
		}

		// if european option
		if (op.type == European) {
			for (int k = N - 2;k >= 0;k--) {
				for (int i = 0;i <= 2 * k;i++) {
					V[i] = disc * (pu*V[i] + pm * V[i + 1] + pd * V[i + 2]);
					if (k <= 2)store[k][i] = V[i];
				}
			}
		}

		// if american option
		else if (op.type == American) {
			for (int k = N - 2;k >= 0;k--) {
				double Smax2 = op.S0*std::pow(u, k);
				for (int i = 0;i <= 2 * k;i++) {
					V[i] = std::max(disc * (pu*V[i] + pm * V[i + 1] + pd * V[i + 2]),op.payoff(Smax2));
					
					if (k <= 2)store[k][i] = V[i];
					Smax2 /= u;
				}
			}
		}
		
		Treevalue=V[0];
		Treedelta = (store[1][0] - store[1][2]) / (op.S0*u - op.S0*d);
		Treegamma = ((store[2][0] - store[2][2]) / (op.S0*u*u - op.S0)
					- (store[2][2] - store[2][4]) / (op.S0 - op.S0*d*d)) / (op.S0*u - op.S0*d);
		Treetheta = (store[1][1] - store[0][0]) /  dt;
		return;
	}
	
	void TBSR(){
		TrinomialTreePricer a(op,N);
		TrinomialTreePricer b(op,N/2);
		a.BSTrinomial();
		b.BSTrinomial();
		Treevalue=(2*a.Treevalue-b.Treevalue);
		Treedelta=(2*a.Treedelta-b.Treedelta);
		Treegamma=(2*a.Treegamma-b.Treegamma);
		Treetheta=(2*a.Treetheta-b.Treetheta);
		return;
	}
	
};