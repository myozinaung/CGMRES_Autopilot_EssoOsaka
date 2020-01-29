#include "nmpc_model.hpp"
// For Single C/GMRES
#include "continuation_gmres.hpp"
#include "cgmres_simulator.hpp"
// For Multi C/GMRES
#include "multiple_shooting_cgmres.hpp"
#include "multiple_shooting_cgmres_simulator.hpp"
// For Multi C/GMRES with Input Saturation
#include "control_input_saturation_sequence.hpp"
#include "multiple_shooting_cgmres_with_saturation.hpp"
#include "multiple_shooting_cgmres_with_saturation_simulator.hpp"

int main()
{
	// Define Simulation Parameters
	double tF  = 200;
	double dt  = 0.1;
	int solver = 1; // 1 = Single C/GMRES, 2 = Multi C/GMRES 
	//(solver 3 needs to change hamiltonian eqs and dimensions)
	std::string filename = "esso_osaka";

	// Initial Conditions
	double x_pos0 = 0;
	double y_pos0 = 0;
	double psi0   = 0 * (pi / 180);
	double u_vel0 = 0.4;
	double v_vel0 = 0;
	double r0     = 0;

	double delta_guess  = 0.01;

	double T_f          = 15; // Prediction Horizon (sec)
	double alpha        = 1.0;
	int horizon_divs    = 20;
	double FD_step      = 1e-6;
	double zeta         = 10; // = (1 / dt)
	int kmax            = 10; // GMRES max iteration no.

	// Set the initial state.
	double initial_state[6] = { x_pos0, y_pos0, psi0, u_vel0, v_vel0, r0 };

	// Set the initial guess of the solution.
	double initial_guess_solution[3] = { delta_guess, delta_guess, delta_guess };


	// Define the model in NMPC.
    NMPCModel nmpc_model;
    // Define the solver.
	
	ContinuationGMRES nmpc_solver(T_f, alpha, horizon_divs, FD_step, zeta, kmax);
	// Initialize the solution of the C/GMRES method.
	nmpc_solver.initSolution(0, initial_state, initial_guess_solution, 1e-06, 100); // tol_res, max_itr (for Newton iteration)

	if (solver == 2) {
		MultipleShootingCGMRES nmpc_solver(T_f, alpha, horizon_divs, FD_step, zeta, kmax);
		std::cout << "Multiple Shooting C/GMRES Solver selected." << std::endl;
		filename = "esso_osaka_multi";
	}
	else if (solver == 3) {
		ControlInputSaturationSequence control_input_saturation_seq;
		control_input_saturation_seq.appendControlInputSaturation(0, -35*pi/180, 35*pi/180, 0.1, 0);
		MultipleShootingCGMRESWithSaturation nmpc_solver(control_input_saturation_seq, 1.0, 1.0, 50, 1e-08, 1000, 3);
		// Set the initial guess of the lagrange multiplier for the condensed constraints with respect to the saturation on the function of the control input .
		double initial_guess_solution2[1] = { delta_guess };
		double initial_guess_lagrange_multiplier[1] = { 0.01 };
		nmpc_solver.initSolution(0, initial_state, initial_guess_solution2, initial_guess_lagrange_multiplier, 1e-06, 100);
		filename = "esso_osaka_multi_sat";
	}
	else {
		std::cout << "Single Shooting C/GMRES Solver selected." << std::endl;
		filename = "esso_osaka_single";
	}

    // Perform a numerical simulation.
    nmpcsim::simulation(nmpc_solver, initial_state, 0, tF, dt, filename); // time step cannot be larger then 0.003
    // t0, tf, dt

    return 0;
}
