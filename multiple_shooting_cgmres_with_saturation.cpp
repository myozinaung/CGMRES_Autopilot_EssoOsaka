#include "multiple_shooting_cgmres_with_saturation.hpp"


inline void MultipleShootingCGMRESWithSaturation::addHamiltonianDerivativeWithControlInput(const double* control_input_and_constraints_vec, const double* saturation_lagrange_multiplier_vec, double* optimality_for_control_input_and_constraints_vec)
{
    for(int i=0; i<dim_saturation_; i++){
        optimality_for_control_input_and_constraints_vec[saturation_seq_.index(i)] += (2*control_input_and_constraints_vec[saturation_seq_.index(i)] - saturation_seq_.min(i) - saturation_seq_.max(i)) * saturation_lagrange_multiplier_vec[i];
    }
}


inline void MultipleShootingCGMRESWithSaturation::computeDummyOptimality(const double* dummy_input_vec, const double* saturation_lagrange_multiplier_vec, double* optimality_for_dummy)
{
    for(int i=0; i<dim_saturation_; i++){
        optimality_for_dummy[i] = 2 * (saturation_seq_.quadratic_weight(i) + saturation_lagrange_multiplier_vec[i]) * dummy_input_vec[i] - saturation_seq_.dummy_weight(i);
    }
}


inline void MultipleShootingCGMRESWithSaturation::computeSaturationOptimality(const double* control_input_and_constraint_vec, const double* dummy_input_vec, double* optimality_for_saturation)
{
    for(int i=0; i<dim_saturation_; i++){
        optimality_for_saturation[i] = control_input_and_constraint_vec[saturation_seq_.index(i)] * (control_input_and_constraint_vec[saturation_seq_.index(i)] - saturation_seq_.min(i) - saturation_seq_.max(i)) + saturation_seq_.min(i) * saturation_seq_.max(i) + dummy_input_vec[i] * dummy_input_vec[i];
    }
}


inline void MultipleShootingCGMRESWithSaturation::computeOptimalityErrorforControlInputAndConstraints(const double time_param, const double* state_vec, const double* control_input_and_constraints_seq, double const* const* state_mat, double const* const* lambda_mat, double const* const* saturation_lagrange_multiplier_mat, double* optimality_for_control_input_and_constraints)
{
    // Set and discretize the horizon.
    double horizon_length = horizon_max_length_ * (1.0 - std::exp(- alpha_ * (time_param - initial_time_)));
    double delta_tau = horizon_length / horizon_division_num_;

    // Compute optimality error for contol input and constraints.
    model_.huFunc(time_param, state_vec, control_input_and_constraints_seq, lambda_mat[0], optimality_for_control_input_and_constraints);
    addHamiltonianDerivativeWithControlInput(control_input_and_constraints_seq, saturation_lagrange_multiplier_mat[0], optimality_for_control_input_and_constraints);

    double tau = time_param+delta_tau;
    for(int i=1; i<horizon_division_num_; i++, tau+=delta_tau){
        model_.huFunc(tau, state_mat[i-1], &(control_input_and_constraints_seq[i*dim_control_input_and_constraints_]), lambda_mat[i], &(optimality_for_control_input_and_constraints[i*dim_control_input_and_constraints_]));
        addHamiltonianDerivativeWithControlInput(&(control_input_and_constraints_seq[i*dim_control_input_and_constraints_]), saturation_lagrange_multiplier_mat[i], &(optimality_for_control_input_and_constraints[i*dim_control_input_and_constraints_]));
    }
}


inline void MultipleShootingCGMRESWithSaturation::computeOptimalityErrorforStateAndLambda(const double time_param, const double* state_vec, const double* control_input_and_constraints_seq, double const* const* state_mat, double const* const* lambda_mat, double** optimality_for_state, double** optimality_for_lambda)
{
    // Set and discretize the horizon.
    double horizon_length = horizon_max_length_ * (1.0 - std::exp(- alpha_ * (time_param - initial_time_)));
    double delta_tau = horizon_length / horizon_division_num_;

    // Compute optimality error for state.
    model_.stateFunc(time_param, state_vec, control_input_and_constraints_seq, dx_vec_);
    for(int j=0; j<dim_state_; j++){
        optimality_for_state[0][j] = state_mat[0][j] - state_vec[j] - delta_tau*dx_vec_[j];
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
 


inline void MultipleShootingCGMRESWithSaturation::computeStateAndLambda(const double time_param, const double* state_vec, const double* control_input_and_constraints_seq, double const* const* optimality_for_state, double const* const* optimality_for_lambda, double** state_mat, double** lambda_mat)
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


inline void MultipleShootingCGMRESWithSaturation::computeOptimalityErrorforSaturation(const double* control_input_and_constraints_seq, double const* const* dummy_input_seq, double const* const* saturation_lagrange_multiplier_seq, double** optimality_for_dummy, double** optimality_for_saturation)
{
    for(int i=0; i<horizon_division_num_; i++){
        computeDummyOptimality(dummy_input_seq[i], saturation_lagrange_multiplier_seq[i], optimality_for_dummy[i]);
    }
    for(int i=0; i<horizon_division_num_; i++){
        computeSaturationOptimality(&(control_input_and_constraints_seq[i*dim_control_input_and_constraints_]), dummy_input_seq[i], optimality_for_saturation[i]);
    }
}


inline void MultipleShootingCGMRESWithSaturation::multiplySaturationErrorInverse(const double* control_input_and_constraints_seq, double const* const* dummy_input_seq, double const* const* saturation_lagrange_multiplier_seq, double const* const* multiplied_dummy_input_seq, double const* const* multiplied_lagrange_multiplier_seq, double** resulted_dummy_input_seq, double** resulted_lagrange_multiplier_seq)
{
    for(int i=0; i<horizon_division_num_; i++){
        for(int j=0; j<dim_saturation_; j++){
            resulted_dummy_input_seq[i][j] = multiplied_lagrange_multiplier_seq[i][j]/(2*dummy_input_seq[i][j]);
        }
    }
    for(int i=0; i<horizon_division_num_; i++){
        for(int j=0; j<dim_saturation_; j++){
            resulted_lagrange_multiplier_seq[i][j] = multiplied_dummy_input_seq[i][j]/(2*dummy_input_seq[i][j]) - ((saturation_lagrange_multiplier_seq[i][j] + saturation_seq_.quadratic_weight(j)) * resulted_dummy_input_seq[i][j])/dummy_input_seq[i][j];
        }
    }
}


inline void MultipleShootingCGMRESWithSaturation::computeDummyOptimalityDifference(const double* control_input_and_constraints_seq, double const* const* dummy_input_seq, const double* control_input_and_constraints_update_seq, double** dummy_difference_seq)
{
    for(int i=0; i<horizon_division_num_; i++){
        for(int j=0; j<dim_saturation_; j++){
            dummy_difference_seq[i][j] = ((2*control_input_and_constraints_seq[i*dim_control_input_and_constraints_+saturation_seq_.index(j)] - saturation_seq_.min(j) - saturation_seq_.max(j)) * control_input_and_constraints_update_seq[i*dim_control_input_and_constraints_+saturation_seq_.index(j)]) / (2*dummy_input_seq[i][j]);
        }
    }
}


inline void MultipleShootingCGMRESWithSaturation::computeSaturationOptimalityDifference(const double* control_input_and_constraints_seq, double const* const* dummy_input_seq, double const* const* saturation_lagrange_multiplier_seq, const double* control_input_and_constraints_update_seq, double** saturation_difference_seq)
{
    for(int i=0; i<horizon_division_num_; i++){
        for(int j=0; j<dim_saturation_; j++){
            saturation_difference_seq[i][j] = - ((saturation_lagrange_multiplier_seq[i][j] + saturation_seq_.quadratic_weight(j)) * (2*control_input_and_constraints_seq[i*dim_control_input_and_constraints_+saturation_seq_.index(j)] - saturation_seq_.min(j) - saturation_seq_.max(j)) * control_input_and_constraints_update_seq[i*dim_control_input_and_constraints_+saturation_seq_.index(j)]) / (2 * dummy_input_seq[i][j] * dummy_input_seq[i][j]);
        }
    }
}


void MultipleShootingCGMRESWithSaturation::bFunc(const double time_param, const double* state_vec, const double* current_solution_vec, double* b_vec)
{
    computeOptimalityErrorforControlInputAndConstraints(time_param, state_vec, current_solution_vec, state_mat_, lambda_mat_, saturation_lagrange_multiplier_mat_, control_input_and_constraints_error_seq_);
    computeOptimalityErrorforControlInputAndConstraints(incremented_time_, incremented_state_vec_, current_solution_vec, state_mat_, lambda_mat_, saturation_lagrange_multiplier_mat_, control_input_and_constraints_error_seq_1_);

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
    computeOptimalityErrorforStateAndLambda(incremented_time_, incremented_state_vec_, current_solution_vec, state_mat_, lambda_mat_, state_error_mat_1_, lambda_error_mat_1_);

    computeOptimalityErrorforSaturation(current_solution_vec, dummy_input_mat_, saturation_lagrange_multiplier_mat_, dummy_error_mat_, saturation_error_mat_);
    for(int i=0; i<horizon_division_num_; i++){
        for(int j=0; j<dim_saturation_; j++){
            dummy_input_mat_1_[i][j] = -zeta_*dummy_error_mat_[i][j];
        }
    }
    for(int i=0; i<horizon_division_num_; i++){
        for(int j=0; j<dim_saturation_; j++){
            saturation_lagrange_multiplier_mat_1_[i][j] = -zeta_*saturation_error_mat_[i][j];
        }
    }
    multiplySaturationErrorInverse(current_solution_vec, dummy_input_mat_, saturation_lagrange_multiplier_mat_, dummy_input_mat_1_, saturation_lagrange_multiplier_mat_1_, dummy_error_mat_1_, saturation_error_mat_1_);
    for(int i=0; i<horizon_division_num_; i++){
        for(int j=0; j<dim_saturation_; j++){
            incremented_saturation_lagrange_multiplier_mat_[i][j] = saturation_lagrange_multiplier_mat_[i][j] + difference_increment_*saturation_error_mat_1_[i][j];
        }
    }
    computeOptimalityErrorforControlInputAndConstraints(incremented_time_, incremented_state_vec_, current_solution_vec, incremented_state_mat_, incremented_lambda_mat_, incremented_saturation_lagrange_multiplier_mat_, control_input_and_constraints_error_seq_3_);

    for(int i=0; i<dim_control_input_and_constraints_seq_; i++){
        incremented_control_input_and_constraints_seq_[i] = current_solution_vec[i] + difference_increment_*control_input_and_constraints_update_seq_[i];
    }
    computeStateAndLambda(incremented_time_, incremented_state_vec_, incremented_control_input_and_constraints_seq_, state_error_mat_1_, lambda_error_mat_1_, incremented_state_mat_, incremented_lambda_mat_);
    computeSaturationOptimalityDifference(current_solution_vec, dummy_input_mat_, saturation_lagrange_multiplier_mat_, control_input_and_constraints_update_seq_, saturation_update_mat_);
    for(int i=0; i<horizon_division_num_; i++){
        for(int j=0; j<dim_saturation_; j++){
            incremented_saturation_lagrange_multiplier_mat_[i][j] = saturation_lagrange_multiplier_mat_[i][j] - difference_increment_*saturation_update_mat_[i][j];
        }
    }
    computeOptimalityErrorforControlInputAndConstraints(incremented_time_, incremented_state_vec_, incremented_control_input_and_constraints_seq_, incremented_state_mat_, incremented_lambda_mat_, incremented_saturation_lagrange_multiplier_mat_, control_input_and_constraints_error_seq_2_);

    for(int i=0; i<dim_control_input_and_constraints_seq_; i++){
        b_vec[i] = (1/difference_increment_ - zeta_) * control_input_and_constraints_error_seq_[i] - control_input_and_constraints_error_seq_3_[i]/difference_increment_ - (control_input_and_constraints_error_seq_2_[i]-control_input_and_constraints_error_seq_1_[i])/difference_increment_;
    }
}


void MultipleShootingCGMRESWithSaturation::axFunc(const double time_param, const double* state_vec, const double* current_solution_vec, const double* direction_vec, double* ax_vec)
{

    for(int i=0; i<dim_control_input_and_constraints_seq_; i++){
        incremented_control_input_and_constraints_seq_[i] = current_solution_vec[i] + difference_increment_*direction_vec[i];
    }
    computeStateAndLambda(incremented_time_, incremented_state_vec_, incremented_control_input_and_constraints_seq_, state_error_mat_1_, lambda_error_mat_1_, incremented_state_mat_, incremented_lambda_mat_);
    computeSaturationOptimalityDifference(current_solution_vec, dummy_input_mat_, saturation_lagrange_multiplier_mat_, direction_vec, saturation_update_mat_);
    for(int i=0; i<horizon_division_num_; i++){
        for(int j=0; j<dim_saturation_; j++){
            incremented_saturation_lagrange_multiplier_mat_[i][j] = saturation_lagrange_multiplier_mat_[i][j] - difference_increment_*saturation_update_mat_[i][j];
        }
    }
    computeOptimalityErrorforControlInputAndConstraints(incremented_time_, incremented_state_vec_, incremented_control_input_and_constraints_seq_, incremented_state_mat_, incremented_lambda_mat_, incremented_saturation_lagrange_multiplier_mat_, control_input_and_constraints_error_seq_2_);

    for(int i=0; i<dim_control_input_and_constraints_seq_; i++){
        ax_vec[i] = (control_input_and_constraints_error_seq_2_[i] - control_input_and_constraints_error_seq_1_[i])/difference_increment_;
    }
}


MultipleShootingCGMRESWithSaturation::MultipleShootingCGMRESWithSaturation(const ControlInputSaturationSequence saturation_seq, const double horizon_max_length, const double alpha, const int horizon_division_num, const double difference_increment, const double zeta, const int dim_krylov) : MatrixFreeGMRES(), 
    model_(), 
    saturation_seq_(saturation_seq), 
    dim_state_(model_.dimState()), 
    dim_control_input_(model_.dimControlInput()), 
    dim_constraints_(model_.dimConstraints()), 
    dim_control_input_and_constraints_(model_.dimControlInput()+model_.dimConstraints()), 
    dim_control_input_and_constraints_seq_(horizon_division_num*(model_.dimControlInput()+model_.dimConstraints())), 
    dim_saturation_(saturation_seq.dimSaturation()), 
    dim_saturation_seq_(horizon_division_num*saturation_seq.dimSaturation()), 
    horizon_division_num_(horizon_division_num), 
    dim_krylov_(dim_krylov), 
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
    state_mat_1_(linearfunc::newMatrix(horizon_division_num, dim_state_)), 
    lambda_mat_(linearfunc::newMatrix(horizon_division_num, dim_state_)), 
    lambda_mat_1_(linearfunc::newMatrix(horizon_division_num, dim_state_)), 
    incremented_state_mat_(linearfunc::newMatrix(horizon_division_num, dim_state_)), 
    incremented_lambda_mat_(linearfunc::newMatrix(horizon_division_num, dim_state_)), 
    state_error_mat_(linearfunc::newMatrix(horizon_division_num, dim_state_)), 
    state_error_mat_1_(linearfunc::newMatrix(horizon_division_num, dim_state_)), 
    lambda_error_mat_(linearfunc::newMatrix(horizon_division_num, dim_state_)), 
    lambda_error_mat_1_(linearfunc::newMatrix(horizon_division_num, dim_state_)), 
    dummy_input_mat_(linearfunc::newMatrix(horizon_division_num, dim_saturation_)), 
    dummy_input_mat_1_(linearfunc::newMatrix(horizon_division_num, dim_saturation_)), 
    saturation_lagrange_multiplier_mat_(linearfunc::newMatrix(horizon_division_num, dim_saturation_)), 
    saturation_lagrange_multiplier_mat_1_(linearfunc::newMatrix(horizon_division_num, dim_saturation_)), 
    incremented_saturation_lagrange_multiplier_mat_(linearfunc::newMatrix(horizon_division_num, dim_saturation_)), 
    dummy_error_mat_(linearfunc::newMatrix(horizon_division_num, dim_saturation_)), 
    dummy_error_mat_1_(linearfunc::newMatrix(horizon_division_num, dim_saturation_)), 
    saturation_error_mat_(linearfunc::newMatrix(horizon_division_num, dim_saturation_)), 
    saturation_error_mat_1_(linearfunc::newMatrix(horizon_division_num, dim_saturation_)), 
    dummy_update_mat_(linearfunc::newMatrix(horizon_division_num, dim_saturation_)), 
    saturation_update_mat_(linearfunc::newMatrix(horizon_division_num, dim_saturation_))
{
    // Set dimensions and parameters in GMRES.
    setGMRESParams((model_.dimControlInput()+model_.dimConstraints())*horizon_division_num, dim_krylov);
}


MultipleShootingCGMRESWithSaturation::~MultipleShootingCGMRESWithSaturation()
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
    linearfunc::deleteMatrix(state_mat_1_);
    linearfunc::deleteMatrix(lambda_mat_);
    linearfunc::deleteMatrix(lambda_mat_1_);
    linearfunc::deleteMatrix(incremented_state_mat_);
    linearfunc::deleteMatrix(incremented_lambda_mat_);
    linearfunc::deleteMatrix(state_error_mat_);
    linearfunc::deleteMatrix(state_error_mat_1_);
    linearfunc::deleteMatrix(lambda_error_mat_);
    linearfunc::deleteMatrix(lambda_error_mat_1_);
    linearfunc::deleteMatrix(dummy_input_mat_);
    linearfunc::deleteMatrix(dummy_input_mat_1_);
    linearfunc::deleteMatrix(saturation_lagrange_multiplier_mat_);
    linearfunc::deleteMatrix(saturation_lagrange_multiplier_mat_1_);
    linearfunc::deleteMatrix(incremented_saturation_lagrange_multiplier_mat_);
    linearfunc::deleteMatrix(dummy_error_mat_);
    linearfunc::deleteMatrix(dummy_error_mat_1_);
    linearfunc::deleteMatrix(saturation_error_mat_);
    linearfunc::deleteMatrix(saturation_error_mat_1_);
    linearfunc::deleteMatrix(dummy_update_mat_);
    linearfunc::deleteMatrix(saturation_update_mat_);
}


void MultipleShootingCGMRESWithSaturation::initSolution(const double initial_time, const double* initial_state_vec, const double* initial_guess_input_vec, const double convergence_radius, const int max_iteration)
{
 //   double initial_solution_vec[dim_control_input_and_constraints_+2*dim_saturation_], initial_lambda_vec[dim_state_], initial_guess_lagrange_multiplier_vec[dim_saturation_];
	double initial_solution_vec[1 + 2 * 1], initial_lambda_vec[4], initial_guess_lagrange_multiplier_vec[1];
    InitCGMRESWithSaturation initializer(saturation_seq_, difference_increment_, dim_krylov_);

    // Intialize the solution
    initial_time_ = initial_time;
    for(int i=0; i<dim_saturation_; i++){
        initial_guess_lagrange_multiplier_vec[i] = 0;
    }
    initializer.solve0stepNOCP(initial_time, initial_state_vec, initial_guess_input_vec, initial_guess_lagrange_multiplier_vec, convergence_radius, max_iteration, initial_solution_vec);

    model_.phixFunc(initial_time, initial_state_vec, initial_lambda_vec);
    for(int i=0; i<horizon_division_num_; i++){
        for(int j=0; j<dim_control_input_and_constraints_; j++){
            control_input_and_constraints_seq_[j] = initial_solution_vec[j];
        }
        for(int j=0; j<dim_saturation_; j++){
            dummy_input_mat_[i][j] = initial_solution_vec[dim_control_input_and_constraints_+j];
        }
        for(int j=0; j<dim_saturation_; j++){
            saturation_lagrange_multiplier_mat_[i][j] = initial_solution_vec[dim_control_input_and_constraints_+dim_saturation_+j];
        }
        for(int j=0; j<dim_state_; j++){
            state_mat_[i][j] = initial_state_vec[j];
        }
        for(int j=0; j<dim_state_; j++){
            lambda_mat_[i][j] = initial_lambda_vec[j];
        }
    }

    // Intialize the optimality error.
 //   double initial_control_input_and_constraints_error[dim_control_input_and_constraints_], initial_dummy_input_error[dim_saturation_], initial_saturation_error[dim_saturation_];
	double initial_control_input_and_constraints_error[1], initial_dummy_input_error[1], initial_saturation_error[1];

    initializer.getControlInputAndConstraintsError(initial_time, initial_state_vec, initial_solution_vec, initial_control_input_and_constraints_error);
    initializer.getDummyInputError(initial_time, initial_state_vec, initial_solution_vec, initial_dummy_input_error);
    initializer.getControlInputSaturationError(initial_time, initial_state_vec, initial_solution_vec, initial_saturation_error);

    for(int i=0; i<horizon_division_num_; i++){
        for(int j=0; j<dim_control_input_and_constraints_; j++){
            control_input_and_constraints_error_seq_[i*dim_control_input_and_constraints_+j] = initial_control_input_and_constraints_error[j];
        }
        for(int j=0; j<dim_saturation_; j++){
            dummy_error_mat_[i][j] = initial_dummy_input_error[j];
        }
        for(int j=0; j<dim_saturation_; j++){
            saturation_error_mat_[i][j] = initial_saturation_error[j];
        }
    }
}


void MultipleShootingCGMRESWithSaturation::initSolution(const double initial_time, const double* initial_state_vec, const double* initial_guess_input_vec, const double initial_guess_lagrange_multiplier, const double convergence_radius, const int max_iteration)
{
 //   double initial_solution_vec[dim_control_input_and_constraints_+2*dim_saturation_], initial_lambda_vec[dim_state_], initial_guess_lagrange_multiplier_vec[dim_saturation_];
	double initial_solution_vec[1 + 2 * 1], initial_lambda_vec[4], initial_guess_lagrange_multiplier_vec[1];
    InitCGMRESWithSaturation initializer(saturation_seq_, difference_increment_, dim_krylov_);

    // Intialize the solution
    initial_time_ = initial_time;
    for(int i=0; i<dim_saturation_; i++){
        initial_guess_lagrange_multiplier_vec[i] = initial_guess_lagrange_multiplier;
    }
    initializer.solve0stepNOCP(initial_time, initial_state_vec, initial_guess_input_vec, initial_guess_lagrange_multiplier_vec, convergence_radius, max_iteration, initial_solution_vec);

    model_.phixFunc(initial_time, initial_state_vec, initial_lambda_vec);
    for(int i=0; i<horizon_division_num_; i++){
        for(int j=0; j<dim_control_input_and_constraints_; j++){
            control_input_and_constraints_seq_[j] = initial_solution_vec[j];
        }
        for(int j=0; j<dim_saturation_; j++){
            dummy_input_mat_[i][j] = initial_solution_vec[dim_control_input_and_constraints_+j];
        }
        for(int j=0; j<dim_saturation_; j++){
            saturation_lagrange_multiplier_mat_[i][j] = initial_solution_vec[dim_control_input_and_constraints_+dim_saturation_+j];
        }
        for(int j=0; j<dim_state_; j++){
            state_mat_[i][j] = initial_state_vec[j];
        }
        for(int j=0; j<dim_state_; j++){
            lambda_mat_[i][j] = initial_lambda_vec[j];
        }
    }

    // Intialize the optimality error.
//    double initial_control_input_and_constraints_error[dim_control_input_and_constraints_], initial_dummy_input_error[dim_saturation_], initial_saturation_error[dim_saturation_];
	double initial_control_input_and_constraints_error[1], initial_dummy_input_error[1], initial_saturation_error[1];

    initializer.getControlInputAndConstraintsError(initial_time, initial_state_vec, initial_solution_vec, initial_control_input_and_constraints_error);
    initializer.getDummyInputError(initial_time, initial_state_vec, initial_solution_vec, initial_dummy_input_error);
    initializer.getControlInputSaturationError(initial_time, initial_state_vec, initial_solution_vec, initial_saturation_error);

    for(int i=0; i<horizon_division_num_; i++){
        for(int j=0; j<dim_control_input_and_constraints_; j++){
            control_input_and_constraints_error_seq_[i*dim_control_input_and_constraints_+j] = initial_control_input_and_constraints_error[j];
        }
        for(int j=0; j<dim_saturation_; j++){
            dummy_error_mat_[i][j] = initial_dummy_input_error[j];
        }
        for(int j=0; j<dim_saturation_; j++){
            saturation_error_mat_[i][j] = initial_saturation_error[j];
        }
    }
}


void MultipleShootingCGMRESWithSaturation::initSolution(const double initial_time, const double* initial_state_vec, const double* initial_guess_input_vec, const double* initial_guess_lagrange_multiplier_vec, const double convergence_radius, const int max_iteration)
{
 //   double initial_solution_vec[dim_control_input_and_constraints_+2*dim_saturation_], initial_lambda_vec[dim_state_];
	double initial_solution_vec[1 + 2 * 1], initial_lambda_vec[4];
    InitCGMRESWithSaturation initializer(saturation_seq_, difference_increment_, dim_krylov_);

    // Intialize the solution
    initial_time_ = initial_time;
    initializer.solve0stepNOCP(initial_time, initial_state_vec, initial_guess_input_vec, initial_guess_lagrange_multiplier_vec, convergence_radius, max_iteration, initial_solution_vec);

    model_.phixFunc(initial_time, initial_state_vec, initial_lambda_vec);
    for(int i=0; i<horizon_division_num_; i++){
        for(int j=0; j<dim_control_input_and_constraints_; j++){
            control_input_and_constraints_seq_[j] = initial_solution_vec[j];
        }
        for(int j=0; j<dim_saturation_; j++){
            dummy_input_mat_[i][j] = initial_solution_vec[dim_control_input_and_constraints_+j];
        }
        for(int j=0; j<dim_saturation_; j++){
            saturation_lagrange_multiplier_mat_[i][j] = initial_solution_vec[dim_control_input_and_constraints_+dim_saturation_+j];
        }
        for(int j=0; j<dim_state_; j++){
            state_mat_[i][j] = initial_state_vec[j];
        }
        for(int j=0; j<dim_state_; j++){
            lambda_mat_[i][j] = initial_lambda_vec[j];
        }
    }

    // Intialize the optimality error.
 //   double initial_control_input_and_constraints_error[dim_control_input_and_constraints_], initial_dummy_input_error[dim_saturation_], initial_saturation_error[dim_saturation_];
	double initial_control_input_and_constraints_error[1], initial_dummy_input_error[1], initial_saturation_error[1];

    initializer.getControlInputAndConstraintsError(initial_time, initial_state_vec, initial_solution_vec, initial_control_input_and_constraints_error);
    initializer.getDummyInputError(initial_time, initial_state_vec, initial_solution_vec, initial_dummy_input_error);
    initializer.getControlInputSaturationError(initial_time, initial_state_vec, initial_solution_vec, initial_saturation_error);

    for(int i=0; i<horizon_division_num_; i++){
        for(int j=0; j<dim_control_input_and_constraints_; j++){
            control_input_and_constraints_error_seq_[i*dim_control_input_and_constraints_+j] = initial_control_input_and_constraints_error[j];
        }
        for(int j=0; j<dim_saturation_; j++){
            dummy_error_mat_[i][j] = initial_dummy_input_error[j];
        }
        for(int j=0; j<dim_saturation_; j++){
            saturation_error_mat_[i][j] = initial_saturation_error[j];
        }
    }
}


void MultipleShootingCGMRESWithSaturation::controlUpdate(const double current_time, const double sampling_period, const double* current_state_vec, double* optimal_control_input_vec)
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
            state_mat_[i][j] += sampling_period * (incremented_state_mat_[i][j] - state_mat_[i][j])/difference_increment_;;
        }
    }
    // lambda_mat_ += sampling_period * (incremented_lambda_mat_ - lambda_mat_)/difference_increment_;
    for(int i=0; i<horizon_division_num_; i++){
        for(int j=0; j<dim_state_; j++){
            lambda_mat_[i][j] += sampling_period * (incremented_lambda_mat_[i][j] - lambda_mat_[i][j])/difference_increment_;
        }
    }

    // Update dummy_input_mat_ and saturation_lagrange_multiplier_mat_.
    computeDummyOptimalityDifference(control_input_and_constraints_seq_, dummy_input_mat_, control_input_and_constraints_update_seq_, dummy_update_mat_);
    computeSaturationOptimalityDifference(control_input_and_constraints_seq_, dummy_input_mat_, saturation_lagrange_multiplier_mat_, control_input_and_constraints_update_seq_, saturation_update_mat_);
    for(int i=0; i<horizon_division_num_; i++){
        for(int j=0; j<dim_saturation_; j++){
            dummy_input_mat_[i][j] += sampling_period * (dummy_error_mat_1_[i][j] - dummy_update_mat_[i][j]);
        }
    }
    for(int i=0; i<horizon_division_num_; i++){
        for(int j=0; j<dim_saturation_; j++){
            saturation_lagrange_multiplier_mat_[i][j] += sampling_period * (saturation_error_mat_1_[i][j] - saturation_update_mat_[i][j]);
        }
    }

    // Update control_input_and_constraints_seq_.
    for(int i=0; i<dim_control_input_and_constraints_seq_; i++){
        control_input_and_constraints_seq_[i] += sampling_period * control_input_and_constraints_update_seq_[i];
    }
    for(int i=0; i<dim_control_input_; i++){
        optimal_control_input_vec[i] = control_input_and_constraints_seq_[i];
    }
}


void MultipleShootingCGMRESWithSaturation::getControlInput(double* control_input_vec) const
{
    for(int i=0; i<dim_control_input_; i++){
        control_input_vec[i] = control_input_and_constraints_seq_[i];
    }
}


double MultipleShootingCGMRESWithSaturation::getError(const double current_time, const double* current_state_vec)
{
    double *control_input_and_constraints_error_seq;
    control_input_and_constraints_error_seq = linearfunc::newVector(dim_control_input_and_constraints_seq_);
    computeOptimalityErrorforControlInputAndConstraints(current_time, current_state_vec, control_input_and_constraints_seq_, state_mat_, lambda_mat_, saturation_lagrange_multiplier_mat_, control_input_and_constraints_error_seq);

    double **state_error_mat, **lambda_error_mat;
    state_error_mat = linearfunc::newMatrix(horizon_division_num_, dim_state_);
    lambda_error_mat = linearfunc::newMatrix(horizon_division_num_, dim_state_);
    computeOptimalityErrorforStateAndLambda(current_time, current_state_vec, control_input_and_constraints_seq_, state_mat_, lambda_mat_, state_error_mat, lambda_error_mat);

    double **dummy_error_mat, **saturation_error_mat;
    dummy_error_mat = linearfunc::newMatrix(horizon_division_num_, dim_saturation_);
    saturation_error_mat = linearfunc::newMatrix(horizon_division_num_, dim_saturation_);
    computeOptimalityErrorforSaturation(control_input_and_constraints_seq_, dummy_input_mat_, saturation_lagrange_multiplier_mat_, dummy_error_mat, saturation_error_mat);

    double squared_error = linearfunc::squaredNorm(dim_control_input_and_constraints_seq_, control_input_and_constraints_error_seq);    
    for(int i=0; i<horizon_division_num_; i++){
        squared_error += (linearfunc::squaredNorm(dim_state_, state_error_mat[i]) + linearfunc::squaredNorm(dim_state_, lambda_error_mat[i]) + linearfunc::squaredNorm(dim_saturation_, dummy_error_mat[i]) + linearfunc::squaredNorm(dim_saturation_, saturation_error_mat[i]));
    }

    linearfunc::deleteVector(control_input_and_constraints_error_seq);
    linearfunc::deleteMatrix(state_error_mat);
    linearfunc::deleteMatrix(lambda_error_mat);
    linearfunc::deleteMatrix(dummy_error_mat);
    linearfunc::deleteMatrix(saturation_error_mat);

    return std::sqrt(squared_error);
}
