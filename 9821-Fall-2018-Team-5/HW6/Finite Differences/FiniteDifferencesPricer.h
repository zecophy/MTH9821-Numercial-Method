#pragma once
#include"Option.h"

// Only for Put!
class FiniteDifferencePircer {
public:
	Option op;
	double x_compute;
	double x_left;
	double x_right;

	double alpha_temp;
	double alpha;
	int M;
	int N;

	double delta_tau;
	double delta_x;
	double tau_final;

	double a;
	double b;

	vector<vector<double>> Domin;// domin[n][m], first is space, second is time 
	vector<double> x_pos;
	vector<double> tau_pos;

	vector<vector<bool>> early_ex;
	
	double V_exact;

	FiniteDifferencePircer(Option op_, double alpha_temp_, int M_, double value_=0.0)
		:op(op_), M(M_), alpha_temp(alpha_temp_) {

		if (op.type == European)
			V_exact = op.BSvalue;
		else if (op.type == American)
			V_exact = value_;

		// Domin parameter
		x_compute = log(op.S0 / op.K);
		tau_final = op.T*op.sigma*op.sigma / 2.0;
		double temp = log(op.S0 / op.K) + (op.r - op.q - op.sigma*op.sigma / 2.0)*op.T;
		x_left = temp - 3.0*op.sigma*sqrt(op.T);
		x_right = temp + 3.0*op.sigma*sqrt(op.T);

		delta_tau = tau_final / (double)M;
		N = int((x_right - x_left) / sqrt(delta_tau / alpha_temp));
		delta_x = (x_right - x_left) / (double)N;
		alpha = delta_tau / delta_x / delta_x;

		// change variable a, b
		a = (op.r - op.q) / op.sigma / op.sigma - 0.5;
		b = ((op.r - op.q) / op.sigma / op.sigma + 0.5)*((op.r - op.q) / op.sigma / op.sigma
			+ 0.5) + 2.0*op.q / op.sigma / op.sigma;

		// Domin
		for (int i = 0;i <= N;++i) {
			vector<double> v(M + 1);
			Domin.emplace_back(v);
			vector<bool> w(M + 1, false);
			early_ex.emplace_back(w);
		}

		double x_temp = x_left;
		for (int i = 0;i <= N;++i) {
			x_pos.emplace_back(x_temp);
			x_temp += delta_x;
		}

		double tau_temp = 0;
		for (int i = 0;i <= M;++i) {
			tau_pos.emplace_back(tau_temp);
			tau_temp += delta_tau;
		}

		// compute
		computeDomin();
		computeVal();

		return;
	}

	// Bonudary functions
	double f(double x) {
		return op.K*exp(a*x)*max(1 - exp(x), 0.0);
	}

	double gleft(double tau) {
		if (op.type == European)
			return op.K*exp(a*x_left + b * tau)*(exp(-2.0*op.r*tau / op.sigma / op.sigma) -
				exp(x_left - 2.0*op.q*tau / op.sigma / op.sigma));
		else if (op.type == American)
			return op.K*exp(a*x_left + b * tau)*(1 - exp(x_left));
	}

	double gright(double tau) {
		return 0.0;
	}


	// main function
	void computeDomin() {
		// initialize
		for (int n = 0;n <= N;++n) {
			Domin[n][0] = f(x_pos[n]);
		}
		for (int m = 1;m <= M;++m) {
			Domin[0][m] = gleft(tau_pos[m]);
			Domin[N][m] = gright(tau_pos[m]);
		}

		// update rest domin
		for (int m = 0;m < M;++m) {
			for (int n = 1;n < N;++n) {
				if (op.type == European) {
					Domin[n][m + 1] = alpha * Domin[n - 1][m] + (1 - 2.0*alpha)*Domin[n][m] 
						+ alpha * Domin[n + 1][m];
				}
				else if (op.type == American) {
					double early_ex_premium = op.K*exp(a*x_pos[n] + b * tau_pos[m+1])
						*max(1 - exp(x_pos[n]), 0.0);

					Domin[n][m + 1] = alpha * Domin[n - 1][m] + (1 - 2.0*alpha)*Domin[n][m]
						+ alpha * Domin[n + 1][m];

					if (early_ex_premium > Domin[n][m + 1]) {
						early_ex[n][m + 1] = true;
						Domin[n][m + 1] = early_ex_premium;
					}
				}
			}
		}
		return;
	}



	double v_approx;
	double v_approx2;

	double error_pointwise;
	double error_pointwise2;

	void computeVal() {
		int i = int((x_compute - x_left) / delta_x);
		double x_i = x_pos[i];
		double x_ip = x_pos[i + 1];
		double V_i = exp(-a * x_i - b * tau_final)*Domin[i][M];
		double V_ip = exp(-a * x_ip - b * tau_final)*Domin[i+1][M];

		double S_i = op.K*exp(x_i);
		double S_ip = op.K*exp(x_ip);

		v_approx = ((S_ip - op.S0)*V_i + (op.S0 - S_i)*V_ip) / (S_ip - S_i);
		error_pointwise = abs(v_approx - V_exact);

		double u_approx = ((x_ip - x_compute)*Domin[i][M] + (x_compute - x_i)*Domin[i + 1][M])
			/ (x_ip - x_i);
		v_approx2 = exp(-a * x_compute - b * tau_final)*u_approx;
		error_pointwise2 = abs(v_approx2 - V_exact);
		return;
	}



	double RMS;

	void computeRMS() {
		double error_total = 0.0;
		int N_RMS = 0;
		for (int i = 0;i <= N;++i) {
			double Sk = op.K*exp(x_pos[i]);
			double Vk = exp(-a * x_pos[i] - b * tau_final)*Domin[i][M];
			double Vexact = Option(op.type, op.cp, Sk, op.K, op.T, op.sigma, op.q, op.r)
				.BSvalue;
			if (Vexact > 0.00001*op.S0) {
				N_RMS += 1;
				error_total += (Vk - Vexact)*(Vk - Vexact) / Vexact / Vexact;
			}
		}
		RMS = sqrt(error_total / (double)N_RMS);
		return;
	}


	double delta_approx;
	double gamma_approx;
	double theta_approx;

	void computeGreeks() {
		int i = int((x_compute - x_left) / delta_x);
		double x_i = x_pos[i];
		double x_ip = x_pos[i + 1];
		double x_im = x_pos[i - 1];
		double x_ipp = x_pos[i + 2];

		double V_i = exp(-a * x_i - b * tau_final)*Domin[i][M];
		double V_ip = exp(-a * x_ip - b * tau_final)*Domin[i + 1][M];
		double V_im = exp(-a * x_im - b * tau_final)*Domin[i - 1][M];
		double V_ipp = exp(-a * x_ipp - b * tau_final)*Domin[i + 2][M];

		double S_i = op.K*exp(x_i);
		double S_ip = op.K*exp(x_ip);
		double S_im = op.K*exp(x_im);
		double S_ipp = op.K*exp(x_ipp);

		delta_approx = (V_ip - V_i) / (S_ip - S_i);
		gamma_approx = ((V_ipp - V_ip) / (S_ipp - S_ip) - (V_i - V_im) / (S_i - S_im))
			/ ((S_ipp + S_ip) / 2.0 - (S_i + S_im) / 2.0);

		double V_i_dt= exp(-a * x_i - b * (tau_final-delta_tau))*Domin[i][M-1];
		double V_ip_dt= exp(-a * x_ip - b * (tau_final - delta_tau))*Domin[i + 1][M-1];

		double v_approx_dt = ((S_ip - op.S0)*V_i_dt + (op.S0 - S_i)*V_ip_dt) / (S_ip - S_i);
		theta_approx = (v_approx_dt - v_approx) / (2.0*delta_tau / op.sigma / op.sigma);
		
		return;
	}

	double v_approx_red;
	double error_pointwise_red;

	void computeRed() {
		if (op.type == American) {
			Option optemp(European, op.cp, op.S0, op.K, op.T, op.sigma, op.q, op.r);
			double v_approx_euro = FiniteDifferencePircer(optemp, alpha_temp, M).v_approx;
			double v_bs = op.BSvalue;
			v_approx_red = v_approx + (v_bs - v_approx_euro);
			error_pointwise_red = abs(v_approx_red - V_exact);
		}
		return;
	}


	vector<double> early_t;
	vector<double> early_s;

	void computeEarly() {
		if (op.type == American) {
			for (int m = 0;m <= M;++m) {
				double t = op.T - 2 * m*delta_tau / op.sigma / op.sigma;
				early_t.emplace_back(t);
				int n;
				for (n = N;n > 1;--n) {
					if (early_ex[n - 1][m])
						break;
				}
				// n == Nopt+1
				double S_Nopt = op.K*exp(x_pos[n - 1]);
				double S_Nopt_p = op.K*exp(x_pos[n]);
				early_s.emplace_back(0.5*(S_Nopt + S_Nopt_p));
			}
		}
		return;
	}

};
