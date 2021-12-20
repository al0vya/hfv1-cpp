// Solver steps
#include "get_nodal_values.h"
#include "get_modal_values.h"
#include "encoding.h"
#include "thresholding.h"
#include "regularisation.h"
#include "extra_significance.h"
#include "decoding.h"
#include "add_ghost_cells.h"
#include "friction_update.h"
#include "get_wet_dry_cells.h"
#include "get_face_values.h"
#include "get_positivity_preserving_nodes.h"
#include "fluxHLL.h"
#include "get_bar_values.h"
#include "fv1_operator.h"
#include "get_dt_CFL.h"

// Aliases
#include "real.h"

// Structures
#include "SimulationParameters.h"
#include "SolverParameters.h"
#include "BoundaryConditions.h"
#include "FlattenedScaleCoeffs.h"
#include "FlattenedDetails.h"
#include "AssembledSolution.h"
#include "Maxes.h"

// Sim/solver setters
#include "set_simulation_parameters.h"
#include "set_solver_parameters.h"
#include "set_boundary_conditions.h"

// Wrapper functions
#include "load_fine_scale_coefficients.h"
#include "get_max_scale_coefficients.h"

// File input/output
#include "write_solution_to_file.h"

int main
(
	int    argc,
	char** argv
)
{
	if (argc < 5)
	{
		printf
		(
			"Please run in the command line as follows:\n\n"
			"HFV1_cpp.exe <TEST_CASE> <NUM_CELLS> <MAX_REF_LVL> <EPSILON>"
		);

		exit(-1);
	}
	
	int steps = 0;

	int  test_case        = strtol(argv[1], nullptr, 10); //set_test_case();
	int  num_cells        = strtol(argv[2], nullptr, 10); //set_num_cells();
	int  refinement_level = strtol(argv[3], nullptr, 10); //set_max_refinement_lvl();
	real epsilon          = strtof(argv[4], nullptr    ); //set_error_threshold_epsilon();

	clock_t start = clock();

	// =========================================================== //
	// INITIALISATION OF VARIABLES AND INSTANTIATION OF STRUCTURES //
	// =========================================================== //

	// Structures
	SimulationParameters sim_params    = set_simulation_parameters(test_case, num_cells);
	SolverParameters     solver_params = set_solver_parameters(epsilon, refinement_level);
	BoundaryConditions   bcs           = set_boundary_conditions(test_case);

	// Variables
	int num_fine_cells = sim_params.cells * (1 << solver_params.L);

	int scale_coeffs_per_cell = (1 << (solver_params.L + 1)) - 1;
	int num_scale_coeffs = scale_coeffs_per_cell * sim_params.cells;

	int details_per_cell = (1 << solver_params.L) - 1;
	int num_details = details_per_cell * sim_params.cells;

	real dx_coarse = (sim_params.xmax - sim_params.xmin) / sim_params.cells;
	real dx_fine = dx_coarse / (1 << solver_params.L);

	bool first_time_step = true;
	real timeNow = 0;
	real dt = C(1e-4);

	NodalValues          nodal_vals(num_fine_cells + 1);
	AssembledSolution    assem_sol(num_fine_cells + 2);
	FaceValues           face_vals(num_fine_cells + 1);
	StarValues           star_vals(num_fine_cells + 1);
	Fluxes               fluxes(num_fine_cells + 1);
	BarValues            bar_vals(num_fine_cells);
	FlattenedScaleCoeffs scale_coeffs(num_scale_coeffs);
	FlattenedDetails     details(num_details);
	Maxes                maxes = { 0, 0, 0 };

	// Arrays
	real* dx_flattened  = new real[num_scale_coeffs];
	real* x_coords      = new real[num_scale_coeffs];
	int*  level_indices = new int[num_scale_coeffs];
	real* x_coarse      = new real[sim_params.cells + 1]();
	real* norm_details  = new real[num_details];
	int*  sig_details   = new int[num_details](); 
	int*  dry_cells     = new int[num_fine_cells + 2];
	real* eta_temp      = new real[num_fine_cells + 2];
	real* delta_west    = new real[num_fine_cells + 1];
	real* delta_east    = new real[num_fine_cells + 1];

	// =========================================================== //

	get_nodal_values
	(
		nodal_vals, 
		sim_params, 
		num_fine_cells, 
		bcs, 
		dx_fine, 
		test_case
	);

	// coarse mesh
	for (int i = 0; i < sim_params.cells + 1; i++) x_coarse[i] = sim_params.xmin + dx_coarse * i;
	
	get_modal_values
	(
		assem_sol, 
		nodal_vals, 
		sim_params
	);

	load_fine_scale_coefficients
	(
		sim_params, 
		solver_params, 
		scale_coeffs_per_cell, 
		assem_sol, 
		scale_coeffs
	);

	while (timeNow < sim_params.simulationTime)
	{
		timeNow += dt;

		if (timeNow - sim_params.simulationTime > 0)
		{
			timeNow -= dt;
			dt = sim_params.simulationTime - timeNow;
			timeNow += dt;
		}

		maxes = get_max_scale_coefficients
		(
			maxes, 
			num_fine_cells, 
			assem_sol, 
			first_time_step, 
			scale_coeffs
		);

		encoding
		(
			sim_params, 
			scale_coeffs_per_cell, 
			details_per_cell,
			num_details,
			solver_params, 
			sig_details, 
			scale_coeffs, 
			details,
			norm_details,
			first_time_step,
			maxes
		);

		thresholding
		(
			sim_params, 
			solver_params, 
			scale_coeffs_per_cell, 
			details_per_cell, 
			norm_details, 
			sig_details
		);

		regularisation
		(
			sim_params, 
			solver_params, 
			details_per_cell, 
			sig_details
		);

		extra_significance
		(
			sim_params, 
			details_per_cell, 
			solver_params, 
			sig_details, 
			norm_details
		);

		decoding
		(
			assem_sol, 
			sim_params, 
			scale_coeffs_per_cell, 
			details_per_cell, 
			x_coords, 
			x_coarse, 
			dx_flattened, 
			dx_coarse, 
			level_indices, 
			sig_details, 
			solver_params, 
			scale_coeffs, 
			details
		);

		add_ghost_cells
		(
			assem_sol,
			bcs
		);

		if (sim_params.manning > 0) friction_update(assem_sol, sim_params, solver_params, dt);

		get_wet_dry_cells
		(
			dry_cells,
			assem_sol,
			solver_params
		);

		get_face_values
		(
			assem_sol,
			face_vals,
			eta_temp
		);

		get_positivity_preserving_nodes
		(
			assem_sol,
			solver_params,
			face_vals,
			star_vals,
			delta_west,
			delta_east
		);

		fluxHLL
		(
			assem_sol,
			solver_params,
			star_vals,
			fluxes
		);

		get_bar_values
		(
			star_vals,
			assem_sol,
			bar_vals
		);

		fv1_operator
		(
			dry_cells,
			fluxes,
			solver_params,
			bar_vals,
			assem_sol,
			dt
		);

		// CFL time step adjustment 
		dt = 1e9;
		real total_mass = 0;

		get_dt_CFL
		(
			solver_params,
			assem_sol,
			dt,
			total_mass
		);

		printf
		(
			"Length: %d, step: %d, time step: %.15f, mass: %.1f, progress: %.17f%%\n",
			assem_sol.length, ++steps, dt, total_mass, timeNow / sim_params.simulationTime * 100
		);
	}

	write_solution_to_file
	(
		sim_params, 
		assem_sol
	);

	// delete buffers
	delete[] x_coarse;
	delete[] dx_flattened;
	delete[] x_coords;
	delete[] level_indices;

	delete[] norm_details;
	delete[] sig_details;

	delete[] dry_cells;
	delete[] eta_temp;
	delete[] delta_west;
	delete[] delta_east;

	clock_t end = clock();

	real end_time = (real)(end - start) / CLOCKS_PER_SEC * C(1000.0);
	printf("Execution time measured using clock(): %f ms\n", end_time);

	return 0;
}