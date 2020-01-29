#include "multiple_shooting_cgmres.hpp"


inline void MultipleShootingCGMRES::computeOptimalityErrorforControlInputAndConstraints(const double time_param, const double* state_vec, const double* control_input_and_constraints_seq, double const* const* state_mat, double const* const* lambda_mat, double* optimality_for_control_input_and_constraints)
{
    // Set and discretize the horizon.
    double horizon_length = horizon_max_length_ * (1.0 - std::exp(- alpha_ * (time_param - initial_time_)));
    double delta_tau = horizon_length / horizon_division_num_;

    // Compute optimality error for control input and constraints.
    model_.huFunc(time_param, state_vec, control_input_and_constraints_seq, lambda_mat[0], optimality_for_control_input_and_constraints);
    double tau = time_param+delta_tau;
    for(int i=1; i<horizon_division_num_; i++, tau+=delta_tau){
        model_.huFunc(tau, state_mat[i-1], &(control_input_and_constraints_seq[i*dim_control_input_and_constraints_]), lambda_mat[i], &(optimality_for_control_input_and_constraints[i*dim_control_input_and_constraints_]));
    }
}


inline void MultipleShootingCGMRES::computeOptimalityErrorforStateAndLambda(const double time_param, const double* state_vec, const double* control_input_and_constraints_seq, double const* const* state_mat, double const* const* lambda_mat, double** optimality_for_state, double** optimality_for_lambda)
{
    // Set and discretize the horizon.
    double horizon_length = horizon_max_length_ * (1.0 - std::exp(- alpha_ * (time_param - initial_time_)));
    double delta_tau = horizon_length / horizon_division_num_;

    // Compute optimality error for state.
    model_.stateFunc(time_param, state_vec, control_input_and_constraints_seq, dx_vec_);
    for(int i=0; i<dim_state_; i++){
        optimality_for_state[0][i] = state_mat[0][i] - state_vec[i] - delta_tau*dx_vec_[i];
    }
    double tau = time_param + delta_tau;
    for(int i=1; i<horizon_division_num_; i++, tau+=delta_tau){
        model_.stateFunc(tau, state_mat[i-1], &(control_input_and_constraints_seq[i*dim_control_input_and_constraints_]), dx_vec_);
        for(int j=0; j<dim_state_; j++){
            optimality_for_state[i][j] = state_mat[i][j] - state_mat[i-1][j] - delta_tau*dx_vec_[j];
        }
    }

    // Compute optimality error for lambda.
    model_.phixFunc(tau, state_mat[horizon_division_num_-1], dx_vec_);
    for(int i=0; i<dim_state_; i++){
        optimality_for_lambda[horizon_division_num_-1][i] = lambda_mat[horizon_division_num_-1][i] - dx_vec_[i];
    }
    for(int i=horizon_division_num_-1; i>=1; i--, tau-=delta_tau){
        model_.hxFunc(tau, state_mat[i-1], &(control_input_and_constraints_seq[i*dim_control_input_and_constraints_]), lambda_mat[i], dx_vec_);
        for(int j=0; j<dim_state_; j++){
            optimality_for_lambda[i-1][j] = lambda_mat[i-1][j] - lambda_mat[i][j] - delta_tau*dx_vec_[j];
        }
    }
}


inline void MultipleShootingCGMRES::computeStateAndLambda(const double time_param, const double* state_vec, const double* control_input_and_constraints_seq, double const* const* optimality_for_state, double const* const* optimality_for_lambda, double** state_mat, double** lambda_mat)
{
    // Set and discretize the horizon.
    double horizon_length = horizon_max_length_ * (1.0 - std::exp(- alpha_ * (time_param - initial_time_)));
    double delta_tau = horizon_length / horizon_division_num_;

    // Compute the sequence of state under the error for state.
    model_.stateFunc(time_param, state_vec, control_input_and_constraints_seq, dx_vec_);
    for(int i=0; i<dim_state_; i++){
        state_mat[0][i] = state_vec[i] + delta_tau*dx_vec_[i] + optimality_for_state[0][i];
    }
    double tau = time_param + delta_tau;
    for(int i=1; i<horizon_division_num_; i++, tau+=delta_tau){
        model_.stateFunc(tau, state_mat[i-1], &(control_input_and_constraints_seq[i*dim_control_input_and_constraints_]), dx_vec_);
        for(int j=0; j<dim_state_; j++){
            state_mat[i][j] = state_mat[i-1][j] + delta_tau*dx_vec_[j] + optimality_for_state[i][j];
        }
    }

    // Compute the sequence of lambda under the error for lambda.
    model_.phixFunc(tau, state_mat[horizon_division_num_-1], dx_vec_);
    for(int i=0; i<dim_state_; i++){
        lambda_mat[horizon_division_num_-1][i] = dx_vec_[i] + optimality_for_lambda[horizon_division_num_-1][i];
    }
    for(int i=horizon_division_num_-1; i>=1; i--, tau-=delta_tau){
        model_.hxFunc(tau, state_mat[i-1], &(control_input_and_constraints_seq[i*dim_control_input_and_constraints_]), lambda_mat[i], dx_vec_);
        for(int j=0; j<dim_state_; j++){
            lambda_mat[i-1][j] = lambda_mat[i][j] + delta_tau*dx_vec_[j] + optimality_for_lambda[i-1][j];
        }
    }
}


void MultipleShootingCGMRES::bFunc(const double time_param, const double* state_vec, const double* current_solution_vec, double* b_vec)
{
    computeOptimalityErrorforControlInputAndConstraints(time_param, state_vec, current_solution_vec, state_mat_, lambda_mat_, control_input_and_constraints_error_seq_);
    computeOptimalityErrorforControlInputAndConstraints(incremented_time_, incremented_state_vec_, current_solution_vec, state_mat_, lambda_mat_, control_input_and_constraints_error_seq_1_);
    computeOptimalityErrorforStateAndLambda(time_param, state_vec, current_solution_vec, state_mat_, lambda_mat_, state_error_mat_, lambda_error_mat_);

    for(int i=0; i<horizon_division_num_; i++){
        for(int j=0; j<dim_state_; j++){
            state_error_mat_1_[i][j] = (1-difference_increment_*zeta_)*state_error_mat_[i][j];
        }
    }
    for(int i=0; i<horizon_division_num_; i++){
        for(int j=0; j<dim_state_; j++){
            lambda_error_mat_1_[i][j] = (1-difference_increment_*zeta_)*lambda_error_mat_[i][j];
        }
    }
    computeStateAndLambda(incremented_time_, incremented_state_vec_, current_solution_vec, state_error_mat_1_, lambda_error_mat_1_, incremented_state_mat_, incremented_lambda_mat_);
    computeOptimalityErrorforControlInputAndConstraints(incremented_time_, incremented_state_vec_, current_solution_vec, incremented_state_mat_, incremented_lambda_mat_, control_input_and_constraints_error_seq_3_);

    computeOptimalityErrorforStateAndLambda(incremented_time_, incremented_state_vec_, current_solution_vec, state_mat_, lambda_mat_, state_error_mat_1_, lambda_error_mat_1_);

    for(int i=0; i<dim_control_input_and_constraints_seq_; i++){
        incremented_control_input_and_constraints_seq_[i] = current_solution_vec[i] + difference_increment_*control_input_and_constraints_update_seq_[i];
    }
    computeStateAndLambda(incremented_time_, incremented_state_vec_, incremented_control_input_and_constraints_seq_, state_error_mat_1_, lambda_error_mat_1_, incremented_state_mat_, incremented_lambda_mat_);    
    computeOptimalityErrorforControlInputAndConstraints(incremented_time_, incremented_state_vec_, incremented_control_input_and_constraints_seq_, incremented_state_mat_, incremented_lambda_mat_, control_input_and_constraints_error_seq_2_);

    for(int i=0; i<dim_control_input_and_constraints_seq_; i++){
        b_vec[i] = (1/difference_increment_ - zeta_)*control_input_and_constraints_error_seq_[i] - control_input_and_constraints_error_seq_3_[i]/difference_increment_ - (control_input_and_constraints_error_seq_2_[i] - control_input_and_constraints_error_seq_1_[i])/difference_increment_;
    }
}


void MultipleShootingCGMRES::axFunc(const double time_param, const double* state_vec, const double* current_solution_vec, const double* direction_vec, double* ax_vec)
{
    for(int i=0; i<dim_control_input_and_constraints_seq_; i++){
        incremented_control_input_and_constraints_seq_[i] = current_solution_vec[i] + difference_increment_*direction_vec[i];
    }
    computeStateAndLambda(incremented_time_, incremented_state_vec_, incremented_control_input_and_constraints_seq_, state_error_mat_1_, lambda_error_mat_1_, incremented_state_mat_, incremented_lambda_mat_);
    computeOptimalityErrorforControlInputAndConstraints(incremented_time_, incremented_state_vec_, incremented_control_input_and_constraints_seq_, incremented_state_mat_, incremented_lambda_mat_, control_input_and_constraints_error_seq_2_);

    for(int i=0; i<dim_control_input_and_constraints_seq_; i++){
        ax_vec[i] = (control_input_and_constraints_error_seq_2_[i] - control_input_and_constraints_error_seq_1_[i])/difference_increment_;
    }
}


MultipleShootingCGMRES::MultipleShootingCGMRES(const double horizon_max_length, const double alpha, const int horizon_division_num, const double difference_increment, const double zeta, const int max_dim_krylov) : MatrixFreeGMRES(), 
    model_(), 
    dim_state_(model_.dimState()), 
    dim_control_input_(model_.dimControlInput()), 
    dim_constraints_(model_.dimConstraints()), 
    dim_control_input_and_constraints_(model_.dimControlInput()+model_.dimConstraints()), 
    dim_state_and_lambda_(2*model_.dimState()), 
    dim_control_input_and_constraints_seq_(horizon_division_num*dim_control_input_and_constraints_), 
    horizon_division_num_(horizon_division_num), 
    max_dim_krylov_(max_dim_krylov), 
    initial_time_(0), 
    horizon_max_length_(horizon_max_length), 
    alpha_(alpha), 
    zeta_(zeta), 
    difference_increment_(difference_increment),
    incremented_time_(0), 
    dx_vec_(linearfunc::newVector(dim_state_)), 
    incremented_state_vec_(linearfunc::newVector(dim_state_)), 
    control_input_and_constraints_seq_(linearfunc::newVector(dim_control_input_and_constraints_seq_)), 
    incremented_control_input_and_constraints_seq_(linearfunc::newVector(dim_control_input_and_constraints_seq_)), 
    control_input_and_constraints_error_seq_(linearfunc::newVector(dim_control_input_and_constraints_seq_)), 
    control_input_and_constraints_error_seq_1_(linearfunc::newVector(dim_control_input_and_constraints_seq_)), 
    control_input_and_constraints_error_seq_2_(linearfunc::newVector(dim_control_input_and_constraints_seq_)), 
    control_input_and_constraints_error_seq_3_(linearfunc::newVector(dim_control_input_and_constraints_seq_)), 
    control_input_and_constraints_update_seq_(linearfunc::newVector(dim_control_input_and_constraints_seq_)), 
    state_mat_(linearfunc::newMatrix(horizon_division_num, dim_state_)), 
    lambda_mat_(linearfunc::newMatrix(horizon_division_num, dim_state_)), 
    incremented_state_mat_(linearfunc::newMatrix(horizon_division_num, dim_state_)), 
    incremented_lambda_mat_(linearfunc::newMatrix(horizon_division_num, dim_state_)), 
    state_error_mat_(linearfunc::newMatrix(horizon_division_num, dim_state_)), 
    state_error_mat_1_(linearfunc::newMatrix(horizon_division_num, dim_state_)), 
    lambda_error_mat_(linearfunc::newMatrix(horizon_division_num, dim_state_)), 
    lambda_error_mat_1_(linearfunc::newMatrix(horizon_division_num, dim_state_))
{
    // Set dimensions and parameters in GMRES.
    setGMRESParams(dim_control_input_and_constraints_seq_, max_dim_krylov);
}


MultipleShootingCGMRES::~MultipleShootingCGMRES()
{
    linearfunc::deleteVector(dx_vec_);
    linearfunc::deleteVector(incremented_state_vec_);
    linearfunc::deleteVector(control_input_and_constraints_seq_);
    linearfunc::deleteVector(incremented_control_input_and_constraints_seq_);
    linearfunc::deleteVector(control_input_and_constraints_error_seq_);
    linearfunc::deleteVector(control_input_and_constraints_error_seq_1_);
    linearfunc::deleteVector(control_input_and_constraints_error_seq_2_);
    linearfunc::deleteVector(control_input_and_constraints_error_seq_3_);
    linearfunc::deleteVector(control_input_and_constraints_update_seq_);
    linearfunc::deleteMatrix(state_mat_);
    linearfunc::deleteMatrix(lambda_mat_);
    linearfunc::deleteMatrix(incremented_state_mat_);
    linearfunc::deleteMatrix(incremented_lambda_mat_);
    linearfunc::deleteMatrix(state_error_mat_);
    linearfunc::deleteMatrix(state_error_mat_1_);
    linearfunc::deleteMatrix(lambda_error_mat_);
    linearfunc::deleteMatrix(lambda_error_mat_1_);
}


void MultipleShootingCGMRES::initSolution(const double initial_time, const double* initial_state_vec, const double* initial_guess_input_vec, const double convergence_radius, const int max_iteration)
{
//    double initial_control_input_and_constraints_vec[dim_control_input_and_constraints_], initial_control_input_and_constraints_error[dim_control_input_and_constraints_], initial_lambda_vec[dim_state_];
	double initial_control_input_and_constraints_vec[6], initial_control_input_and_constraints_error[6], initial_lambda_vec[12];
    InitCGMRES initializer(difference_increment_, max_dim_krylov_);
    initial_time_ = initial_time;

    // Intialize the solution
    initializer.solve0stepNOCP(initial_time, initial_state_vec, initial_guess_input_vec, convergence_radius, max_iteration, initial_control_input_and_constraints_vec);
    model_.phixFunc(initial_time, initial_state_vec, initial_lambda_vec);
    for(int i=0; i<horizon_division_num_; i++){
        for(int j=0; j<dim_control_input_and_constraints_; j++){
            control_input_and_constraints_seq_[i*dim_control_input_and_constraints_+j] = initial_control_input_and_constraints_vec[j];
        }
        for(int j=0; j<dim_state_; j++){
            state_mat_[i][j] = initial_state_vec[j];
        }
        for(int j=0; j<dim_state_; j++){
            lambda_mat_[i][j] = initial_lambda_vec[j];
        }
    }

    // Intialize the optimality error.
    initializer.getOptimalityErrorVec(initial_time, initial_state_vec, initial_control_input_and_constraints_vec, initial_control_input_and_constraints_error);
    for(int i=0; i<horizon_division_num_; i++){
        for(int j=0; j<dim_control_input_and_constraints_; j++){
            control_input_and_constraints_error_seq_[i*dim_control_input_and_constraints_+j] = initial_control_input_and_constraints_error[j];
        }
    }
}


void MultipleShootingCGMRES::controlUpdate(const double current_time, const double sampling_period, const double* current_state_vec, double* optimal_control_input_vec)
{
    // Predict the incremented state.
    incremented_time_ = current_time + difference_increment_;
    model_.stateFunc(current_time, current_state_vec, control_input_and_constraints_seq_, dx_vec_);
    for(int i=0; i<dim_state_; i++){
        incremented_state_vec_[i] = current_state_vec[i] + difference_increment_*dx_vec_[i];
    }

    forwardDifferenceGMRES(current_time, current_state_vec, control_input_and_constraints_seq_, control_input_and_constraints_update_seq_);

    // Update state_mat_ and lamdba_mat_ by the difference approximation.
    for(int i=0; i<dim_control_input_and_constraints_seq_; i++){
        incremented_control_input_and_constraints_seq_[i] = control_input_and_constraints_seq_[i] + difference_increment_*control_input_and_constraints_update_seq_[i];
    }
    for(int i=0; i<horizon_division_num_; i++){
        for(int j=0; j<dim_state_; j++){
            state_error_mat_1_[i][j] = (1-difference_increment_*zeta_)*state_error_mat_[i][j];
        }
    }
    for(int i=0; i<horizon_division_num_; i++){
        for(int j=0; j<dim_state_; j++){
            lambda_error_mat_1_[i][j] = (1-difference_increment_*zeta_)*lambda_error_mat_[i][j];
        }
    }
    computeStateAndLambda(incremented_time_, incremented_state_vec_, incremented_control_input_and_constraints_seq_, state_error_mat_1_, lambda_error_mat_1_, incremented_state_mat_, incremented_lambda_mat_);

    // state_mat_ += sampling_period * (incremented_state_mat_ - state_mat_)/difference_increment_;
    for(int i=0; i<horizon_division_num_; i++){
        for(int j=0; j<dim_state_; j++){
            state_mat_[i][j] += sampling_period * (incremented_state_mat_[i][j] - state_mat_[i][j])/difference_increment_;
        }
    }
    // lambda_mat_ += sampling_period * (incremented_lambda_mat_ - lambda_mat_)/difference_increment_;
    for(int i=0; i<horizon_division_num_; i++){
        for(int j=0; j<dim_state_; j++){
            lambda_mat_[i][j] += sampling_period * (incremented_lambda_mat_[i][j] - lambda_mat_[i][j])/difference_increment_;
        }
    }

    // Update control_input_and_constraints_seq_
    for(int i=0; i<dim_control_input_and_constraints_seq_; i++){
        control_input_and_constraints_seq_[i] += sampling_period * control_input_and_constraints_update_seq_[i];
    }
    for(int i=0; i<dim_control_input_; i++){
        optimal_control_input_vec[i] = control_input_and_constraints_seq_[i];
    }
}


double MultipleShootingCGMRES::getError(const double current_time, const double* current_state_vec)
{
    computeOptimalityErrorforControlInputAndConstraints(current_time, current_state_vec, control_input_and_constraints_seq_, state_mat_, lambda_mat_, control_input_and_constraints_error_seq_);
    computeOptimalityErrorforStateAndLambda(current_time, current_state_vec, control_input_and_constraints_seq_, state_mat_, lambda_mat_, state_error_mat_, lambda_error_mat_);

    double squared_error = linearfunc::squaredNorm(dim_control_input_and_constraints_seq_, control_input_and_constraints_error_seq_);
    for(int i=0; i<horizon_division_num_; i++){
        squared_error += (linearfunc::squaredNorm(dim_state_, state_error_mat_[i]) + linearfunc::squaredNorm(dim_state_, lambda_error_mat_[i]));
    }
    return std::sqrt(squared_error);
}


void MultipleShootingCGMRES::getControlInput(double* control_input_vec) const
{
    for(int i=0; i<dim_control_input_; i++){
        control_input_vec[i] = control_input_and_constraints_seq_[i];
    }
}