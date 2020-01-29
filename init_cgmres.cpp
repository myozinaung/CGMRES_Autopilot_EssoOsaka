#include "init_cgmres.hpp"


inline void InitCGMRES::computeOptimalityErrors(const double time_param, const double* state_vec, const double* current_solution_vec, double* optimality_error_vec)
{
    model_.phixFunc(time_param, state_vec, lambda_vec_);
    model_.huFunc(time_param, state_vec, current_solution_vec, lambda_vec_, optimality_error_vec);
}


void InitCGMRES::bFunc(const double time_param, const double* state_vec, const double* current_solution_vec, double* b_vec)
{
    for(int i=0; i<dim_solution_; i++){
        incremented_solution_vec_[i] = current_solution_vec[i] + finite_diff_step_*solution_update_vec_[i];
    }
    computeOptimalityErrors(time_param, state_vec, current_solution_vec, error_vec_);
    computeOptimalityErrors(time_param, state_vec, incremented_solution_vec_, error_vec_1_);

    for(int i=0; i<dim_solution_; i++){
        b_vec[i] = - error_vec_[i] - (error_vec_1_[i] - error_vec_[i])/finite_diff_step_;
    }
}


void InitCGMRES::axFunc(const double time_param, const double* state_vec, const double* current_solution_vec, const double* direction_vec, double* ax_vec)
{
    for(int i=0; i<dim_solution_; i++){
        incremented_solution_vec_[i] = current_solution_vec[i] + finite_diff_step_*direction_vec[i];
    }
    computeOptimalityErrors(time_param, state_vec, incremented_solution_vec_, error_vec_1_);

    for(int i=0; i<dim_solution_; i++){
        ax_vec[i] = (error_vec_1_[i] - error_vec_[i])/finite_diff_step_;
    }
}


InitCGMRES::InitCGMRES(const double finite_diff_step, const int kmax) : MatrixFreeGMRES(), 
    model_(), 
    dim_solution_(model_.dimControlInput()+model_.dimConstraints()), 
    finite_diff_step_(finite_diff_step), 
    solution_update_vec_(linearfunc::newVector(dim_solution_)), 
    incremented_solution_vec_(linearfunc::newVector(dim_solution_)), 
    lambda_vec_(linearfunc::newVector(model_.dimState())), 
    error_vec_(linearfunc::newVector(dim_solution_)), 
    error_vec_1_(linearfunc::newVector(dim_solution_)), 
    error_vec_2_(linearfunc::newVector(dim_solution_))
{
    // Set dimensions and parameters in GMRES.
    setGMRESParams(dim_solution_, kmax);
}


InitCGMRES::~InitCGMRES()
{
    linearfunc::deleteVector(solution_update_vec_);
    linearfunc::deleteVector(incremented_solution_vec_);
    linearfunc::deleteVector(lambda_vec_);
    linearfunc::deleteVector(error_vec_);
    linearfunc::deleteVector(error_vec_1_);
    linearfunc::deleteVector(error_vec_2_);
}


void InitCGMRES::solve0stepNOCP(const double initial_time, const double* initial_state_vec, const double* initial_guess_vec, const double convergence_radius, const int max_iteration, double* solution_vec)
{
    for(int i=0; i<dim_solution_; i++){
        solution_vec[i] = initial_guess_vec[i];
    }
    computeOptimalityErrors(initial_time, initial_state_vec, solution_vec, error_vec_);
    int j=0;
    while(linearfunc::squaredNorm(dim_solution_, error_vec_) > convergence_radius*convergence_radius && j < max_iteration){
        forwardDifferenceGMRES(initial_time, initial_state_vec, solution_vec, solution_update_vec_);
        for(int i=0; i<dim_solution_; i++){
            solution_vec[i] += solution_update_vec_[i];
        }
        computeOptimalityErrors(initial_time, initial_state_vec, solution_vec, error_vec_);
        j++;
		//double errVal = linearfunc::squaredNorm(dim_solution_, error_vec_);
    }
}


void InitCGMRES::getOptimalityErrorVec(const double initial_time, const double* initial_state_vec, const double* current_solution_vec, double* error_vec)
{
    computeOptimalityErrors(initial_time, initial_state_vec, current_solution_vec, error_vec);
}