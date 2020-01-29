#include "hamiltonian.hpp"

double hamiltonian::hamiltonian_fun(const double t, const double* x, const double* u, const double* lmd)
{
	double f[dim_state_];
	StateODE ode;
	ode.MMG_EO_ODE(t, x, u, f);

	double L = stageCost(t, x, u);
	double h_constr = constr_fun(t, x, u);
	double h_f = 0;
	for (int i = 0; i < dim_state_; i++) {
		h_f += lmd[i] * f[i];
	}
	double h = L + h_f + h_constr;
	return h;
}

double hamiltonian::phi_fun(const double t, const double* x)
{
	double phi = 0;
	for (int i = 0; i < dim_state_; ++i)
	{
		phi += q_terminal[i] * pow((x[i] - x_ref[i]), 2);
	}
	return phi / 2;
}

double hamiltonian::stageCost(const double t, const double* x, const double* u)
{
	double costStateTrack = 0;
	for (int i = 0; i < dim_state_; ++i)
	{
		costStateTrack += q[i] * pow((x[i] - x_ref[i]), 2);
	}
	double costControlTrack = r[0] * u[0] * u[0];
	double costDummy = dummy_wt * u[1];
	double L = costStateTrack / 2 + costControlTrack / 2 - costDummy; //  cost dummy is NEGATIVE
	return L;
}


double hamiltonian::constr_fun(const double t, const double* x, const double* u)
{
	double C = u[0] * u[0] + u[1] * u[1] - pow(((u_max - u_min) / 2), 2);
	double mu_dummy = u[2];
	double h_constr = mu_dummy * C;
	return h_constr;
}

