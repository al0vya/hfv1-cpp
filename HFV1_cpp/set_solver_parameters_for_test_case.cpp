#include "set_solver_parameters_for_test_case.h"
#include "real.h"

SolverParameters set_solver_parameters_for_test_case(int user_input_max_refinement_level, real user_input_epsilon)
{
	SolverParameters solverParameters;

	solverParameters.CFL = C(0.33);
	solverParameters.tolDry = C(1e-4);
	solverParameters.g = C(9.80665);
	solverParameters.L = user_input_max_refinement_level;
	solverParameters.epsilon = user_input_epsilon;

	return solverParameters;
}