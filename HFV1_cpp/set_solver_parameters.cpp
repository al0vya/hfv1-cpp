#include "set_solver_parameters.h"

SolverParameters set_solver_parameters
(
	real& user_input_epsilon, 
	int&  user_input_max_refinement_level
)
{
	SolverParameters solverParameters;

	solverParameters.CFL = C(0.33);
	solverParameters.tol_dry = C(1e-4);
	solverParameters.g = C(9.80665);
	solverParameters.epsilon = user_input_epsilon;
	solverParameters.L = user_input_max_refinement_level;

	return solverParameters;
}