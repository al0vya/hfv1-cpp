#include "set_simulation_parameters.h"

SimulationParameters set_simulation_parameters(int test_case_selection, int num_cells)
{
	SimulationParameters sim_params;

	switch (test_case_selection)
	{
	case 1:
		sim_params.cells = num_cells;
		sim_params.xmin = 0;
		sim_params.xmax = 50;
		sim_params.simulationTime = C(2.5);
		sim_params.manning = C(0.0);
		break;
	case 2:
		sim_params.cells = num_cells;
		sim_params.xmin = 0;
		sim_params.xmax = 50;
		sim_params.simulationTime = C(1.3);
		sim_params.manning = C(0.0);
		break;
	case 3:
		sim_params.cells = num_cells;
		sim_params.xmin = 0;
		sim_params.xmax = 50;
		sim_params.simulationTime = C(1.3);
		sim_params.manning = C(0.02);
		break;
	case 4:
		sim_params.cells = num_cells;
		sim_params.xmin = 0;
		sim_params.xmax = 50;
		sim_params.simulationTime = C(0.5);
		sim_params.manning = C(0.0);
		break;
	case 5:
		sim_params.cells = num_cells;
		sim_params.xmin = 0;
		sim_params.xmax = 50;
		sim_params.simulationTime = C(0.5);
		sim_params.manning = C(0.0);
		break;
	case 6:
		sim_params.cells = num_cells;
		sim_params.xmin = 0;
		sim_params.xmax = 50;
		sim_params.simulationTime = C(10.0);
		sim_params.manning = C(0.02);
		break;
	default:

		break;
	}

	return sim_params;
}