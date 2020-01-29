#ifndef HAMILTONIAN_H
#define HAMILTONIAN_H

#include "MMG_EO.hpp"
#include "linear_funcs.hpp"

class hamiltonian{
private:

    
public:

    static constexpr int dim_state_         = 6;
    static constexpr int dim_control_input_ = 2; // real u and dummy u
    static constexpr int dim_constraints_   = 1; // mu

    double q[6]          = { 0, 0, 1, 0, 0, 0 };
    double r[2]          = { 0, 0 };
    double q_terminal[6] = { 0, 0, 1, 0, 0, 0 };
    //double x_ref[6]      = { 0, 0, 0 * pi / 180, 0.4, 0, 0 };

    double dummy_wt = 0.01;

    // Control Inputs Constraints
    double u_min = -35 * pi / 180;
    double u_max =  35 * pi / 180;

    double hamiltonian_fun(const double t, const double* x, const double* u, const double* lmd);
    double phi_fun(const double t, const double* x);
    double stageCost(const double t, const double* x, const double* u);
    double constr_fun(const double t, const double* x, const double* u);

};

#endif // !HAMILTONIAN_H