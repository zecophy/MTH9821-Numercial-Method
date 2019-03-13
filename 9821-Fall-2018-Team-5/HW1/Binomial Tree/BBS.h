#pragma once
#include"AverageBinomialTree.h"
#include<cmath>

double cdf(double d) {
	return (1 + std::erf(d/std::sqrt(2))) / 2;
}

ValueAndGreeks BS(Option& op) {
	double d1 = (std::log(op.S0 / op.K) + (op.r - op.q + op.sigma * op.sigma / 2)*op.T) / op.sigma / std::sqrt(op.T);
	double d2 = d1 - op.sigma*std::sqrt(op.T);
	double v=0;
	double delta = 0, gamma = 0, theta = 0;
	double Nd1 = 0, Nd2 = 0;
	double pi = 3.1415926;

	if (op.cp == call) {
		Nd1 = cdf(d1);
		Nd2 = cdf(d2);
		v= op.S0*std::exp(-op.q*op.T)*Nd1 - op.K*std::exp(-op.r*op.T)*Nd2;
		delta = std::exp(-op.q*op.T)*Nd1;
		theta -= op.S0*op.sigma*std::exp(-op.q*op.T) / 2 / std::sqrt(2 * pi*op.T)*std::exp(-d1 * d1 / 2);
		theta += op.q*op.S0*std::exp(-op.q*op.T)*Nd1;
		theta -= op.r*op.K*std::exp(-op.r*op.T)*Nd2;
	}
	else if (op.cp == put) {
		Nd1 = cdf(-d1);
		Nd2 = cdf(-d2);
		v= -op.S0*std::exp(-op.q*op.T)*Nd1 + op.K*std::exp(-op.r*op.T)*Nd2;
		delta = -std::exp(-op.q*op.T)*Nd1;
		theta -= op.S0*op.sigma*std::exp(-op.q*op.T) / 2 / std::sqrt(2 * pi*op.T)*std::exp(-d1 * d1 / 2);
		theta -= op.q*op.S0*std::exp(-op.q*op.T)*Nd1;
		theta += op.r*op.K*std::exp(-op.r*op.T)*Nd2;
	}
	gamma = std::exp(-op.q*op.T) / (op.S0*op.sigma*std::sqrt(op.T)) / std::sqrt(2 * pi)
		*std::exp(-d1 * d1 / 2);

	return { v,delta,gamma,theta };
}

ValueAndGreeks BBS(Option& op, int N = 100) {
	double dt = op.T / N;
	double disc = std::exp(-op.r*dt);
	double u = std::exp(op.sigma*std::sqrt(dt));
	double d = 1 / u;
	double pu = (std::exp((op.r - op.q)*dt) - d) / (u - d);
	double pd = 1 - pu;

	double Smax = op.S0 * std::pow(u, N-1);
	std::vector<double> V(N);

	if (op.type == european)
		for (int i = 0;i < N;i++) {
			Option temp = op;
			temp.T = dt;
			temp.S0 = Smax;
			V[i] = BS(temp).value;
			Smax /= (u*u);
		}

	if (op.type == american)
		for (int i = 0;i < N;i++) {
			Option temp = op;
			temp.T = dt;
			temp.S0 = Smax;
			if (op.cp == call)
				V[i] = std::max(BS(temp).value, -op.K + Smax);
			else if (op.cp == put)
				V[i] = std::max(BS(temp).value, op.K - Smax);
			Smax /= (u*u);
		}

	std::vector<std::vector<double>> store{ { 0 },{ 0,0 },{ 0,0,0 } };

	// if european option
	if (op.type == european) {
		for (int k = N - 2;k >= 0;k--) {
			for (int i = 0;i <= k;i++) {
				V[i] = disc * (pu*V[i] + pd * V[i + 1]);
				if (k <= 2)store[k][i] = V[i];
			}
		}
	}

	// if american option
	else if (op.type == american) {
		for (int k = N - 2;k >= 0;k--) {
			double Smax2 = op.S0*std::pow(u, k);
			for (int i = 0;i <= k;i++) {
				// if put
				if (op.cp == put) {
					V[i] = std::max(disc*(pu*V[i] + pd * V[i + 1]), std::max(op.K - Smax2, 0.0));
				}
				//if call
				else if (op.cp == call) {
					V[i] = std::max(disc*(pu*V[i] + pd * V[i + 1]), std::max(-op.K + Smax2, 0.0));
				}
				if (k <= 2)store[k][i] = V[i];
				Smax2 /= (u*u);
			}
		}
	}
	double delta = (store[1][0] - store[1][1]) / (op.S0*u - op.S0*d);
	double gamma = ((store[2][0] - store[2][1]) / (op.S0*u*u - op.S0)
		- (store[2][1] - store[2][2]) / (op.S0 - op.S0*d*d)) / ((op.S0*u*u - op.S0*d*d) / 2);
	double theta = (store[2][1] - store[0][0]) / 2 / dt;
	return{ V[0],delta,gamma,theta };
}
