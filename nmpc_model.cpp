#include "nmpc_model.hpp"

// State equation f(t, x, u)
// t : time parameter
// x : state vector
// u : control input vector
// f : the value of f(t, x, u)
void NMPCModel::stateFunc(const double t, const double* x, const double* u, double* f)
{
    StateODE ode_fun;
    ode_fun.MMG_EO_ODE(t, x, u, f);
}

// Partial derivative of terminal cost with respect to state, dphi/dx(t, x)
// phi(t, x) = (q_terminal[0]*(x[0]-x_ref[0])^2 + q_terminal[1]*(x[1]-x_ref[1])^2 + q_terminal[2]*(x[2]-x_ref[2])^2 + q_terminal[3]*(x[3]-x_ref[3])^2)/2
// t    : time parameter
// x    : state vector
// phix : the value of dphi/dx(t, x)
void NMPCModel::phixFunc(const double t, const double* x, double* phix) // Analytical Gradient
{
    hamiltonian ham;
    constexpr int dim_state_ = ham.dim_state_;
    for (int i = 0; i < dim_state_; i++) {

		phix[i] = (1.0 / 2.0) * q_terminal[i] * (2 * x[i] - 2 * x_ref[i]);
	}
}

// Partial derivative of the Hamiltonian with respect to state, dH/dx(t, x, u, lmd)
// H(t, x, u, lmd) = L(t, x, u) + lmd * f(t, x, u)
// L(t, x, u) = (q[0]*(x[0]-x_ref[0])^2 + q[1]*(x[1]-x_ref[1])^2 + q[2]*(x[2]-x_ref[2])^2 + q[3]*(x[3]-x_ref[3])^2 + r[0]*u[0]^2)/2
// t   : time parameter
// x   : state vector
// u   : control input vector
// lmd : the Lagrange multiplier for the state equation
// hx  : the value of dH/dx(t, x, u, lmd)
void NMPCModel::hxFunc(const double t, const double* x, const double* u, const double* lmd, double* hx)
{
    hamiltonian ham;
    constexpr int dim_state_ = ham.dim_state_;
    double h1 = 0, h2 = 0, h3 = 0, h4 = 0;
    double x_new1[dim_state_];
    double x_new2[dim_state_];
    double x_new3[dim_state_];
    double x_new4[dim_state_];

    for (int i = 0; i < dim_state_; ++i) {
        for (int j = 0; j < dim_state_; ++j) {
            x_new1[j]  = x[j];
            x_new2[j] = x[j];
            x_new3[j] = x[j];
            x_new4[j] = x[j];
        }

        x_new1[i] = x[i] + dhx;
        x_new2[i] = x[i] - dhx;
        x_new3[i] = x[i] + 2*dhx;
        x_new4[i] = x[i] - 2*dhx;
        h1 = ham.hamiltonian_fun(t, x_new1, u, lmd);
        h2 = ham.hamiltonian_fun(t, x_new2, u, lmd);
        h3 = ham.hamiltonian_fun(t, x_new3, u, lmd);
        h4 = ham.hamiltonian_fun(t, x_new4, u, lmd);
        hx[i] = (h4 - 8*h2 + 8*h1 - h3) / (12 * dhx);
    }
}

// Partial derivative of the Hamiltonian with respect to control input and constraints, dH/du(t, x, u, lmd)
// H(t, x, u, lmd) = L(t, x, u) + lmd * f(t, x, u)
// L(t, x, u) = (q[0]*(x[0]-x_ref[0])^2 + q[1]*(x[1]-x_ref[1])^2 + q[2]*(x[2]-x_ref[2])^2 + q[3]*(x[3]-x_ref[3])^2 + r[0]*u[0]^2)/2
// t   : time parameter
// x   : state vector
// u   : control input vector
// lmd : the Lagrange multiplier for the state equation
// hu  : the value of dH/du(t, x, u, lmd)
void NMPCModel::huFunc(const double t, const double* x, const double* u, const double* lmd, double* hu)
{
    // 3 variables for each control input >> "real u", "dummpy u", "mu"
    hamiltonian ham;
    constexpr int dim_control_input_ = ham.dim_control_input_;
    constexpr int dim_constraints_   = ham.dim_constraints_;
    double h1 = 0, h2 = 0, h3 = 0, h4 = 0;
    double u_new1[dim_control_input_ + dim_constraints_];
    double u_new2[dim_control_input_ + dim_constraints_];
    double u_new3[dim_control_input_ + dim_constraints_];
    double u_new4[dim_control_input_ + dim_constraints_];

    for (int i = 0; i < (dim_control_input_ + dim_constraints_); i++)
    {
        for (int j = 0; j < (dim_control_input_ + dim_constraints_); j++) {
            u_new1[j] = u[j];
            u_new2[j] = u[j];
            u_new3[j] = u[j];
            u_new4[j] = u[j];
        }

        u_new1[i] = u[i] + dhu;
        u_new2[i] = u[i] - dhu;
        u_new3[i] = u[i] + 2*dhu;
        u_new4[i] = u[i] - 2*dhu;
        h1 = ham.hamiltonian_fun(t, x, u_new1, lmd);
        h2 = ham.hamiltonian_fun(t, x, u_new2, lmd);
        h3 = ham.hamiltonian_fun(t, x, u_new3, lmd);
        h4 = ham.hamiltonian_fun(t, x, u_new4, lmd);
        hu[i] = (h4 - 8 * h2 + 8 * h1 - h3) / (12 * dhu);
    }
}
