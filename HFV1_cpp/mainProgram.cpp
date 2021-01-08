#include <stdlib.h>
#include <stdio.h>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <fstream>
#include <time.h>

// Solver steps
#include "get_nodal_values.h"
#include "get_modal_values.h"
#include "add_ghost_cells.h"
#include "friction_update.h"
#include "get_wet_dry_cells.h"
#include "get_face_values.h"
#include "get_positivity_preserving_nodes.h"
#include "fluxHLL.h"
#include "get_bar_values.h"
#include "fv1_operator.h"
#include "get_dt_CFL.h"

#include "real.h"

// Structures
#include "SimulationParameters.h"
#include "SolverParameters.h"
#include "BoundaryConditions.h"
#include "FlattenedScaleCoeffs.h"
#include "FlattenedDetails.h"
#include "AssembledSolution.h"

// Sim/solver setters
#include "set_simulation_parameters.h"
#include "set_solver_parameters.h"
#include "set_boundary_conditions.h"
#include "set_max_refinement_lvl.h"
#include "set_test_case.h"
#include "set_num_cells.h"
#include "set_error_threshold_epsilon.h"

#include "bedDataConservative.h"
#include "hInitialConservative.h"
#include "hInitialOvertopping.h"
#include "qInitial.h"
#include "friction_update.h"
#include "fluxHLL.h"
#include "encodeDetail.h"
#include "encodeScale.h"
#include "decode1.h"
#include "decode2.h"
#include "treeTraversalDecode.h"
#include "write_solution_to_file.h"

// Memory (de)allocators
#include "malloc_fluxes.h"
#include "malloc_assembled_solution.h"
#include "malloc_bar_values.h"
#include "malloc_details.h"
#include "malloc_nodal_values.h"
#include "malloc_star_values.h"
#include "malloc_scale_coefficients.h"
#include "malloc_face_values.h"
#include "free_fluxes.h"
#include "free_assembled_solution.h"
#include "free_bar_values.h"
#include "free_details.h"
#include "free_nodal_values.h"
#include "free_star_values.h"
#include "free_scale_coefficients.h"
#include "free_face_values.h"

int main()
{
	clock_t start = clock();

	// quintessential for-loop index
	int i;
	int steps = 0;

	int  test_case = set_test_case();
	int  num_cells = set_num_cells();
	int  refinement_level = set_max_refinement_lvl();
	real epsilon = set_error_threshold_epsilon();

	// =========================================================== //
	// INITIALISATION OF VARIABLES AND INSTANTIATION OF STRUCTURES //
	// =========================================================== //

	// Structures
	SimulationParameters sim_params    = set_simulation_parameters(test_case, num_cells);
	SolverParameters     solver_params = set_solver_parameters(epsilon, refinement_level);
	BoundaryConditions   bcs           = set_boundary_conditions(test_case);

	NodalValues          nodal_vals;
	AssembledSolution    assem_sol;
	FaceValues           face_vals;
	StarValues           star_vals;
	Fluxes               fluxes;
	BarValues            bar_vals;
	FlattenedScaleCoeffs scale_coeffs;
	FlattenedDetails     details;
	
	// Variables
	int num_fine_cells = sim_params.cells * (1 << solver_params.L);

	int scaleCoeffsPerCell = (1 << (solver_params.L + 1)) - 1;
	int totalScaleCoeffs = scaleCoeffsPerCell * sim_params.cells;

	int detailsPerCell = (1 << solver_params.L) - 1;
	int totalDetails = detailsPerCell * sim_params.cells;

	assem_sol.length = 0;

	real qMax = 0;
	real zMax = 0;
	real etaMax = 0;

	bool firstStep = true;
	real timeNow = 0;
	real dt = C(1e-4);

	// Memory allocation
	malloc_nodal_values(nodal_vals, num_fine_cells);
	malloc_assembled_solution(assem_sol, num_fine_cells);
	malloc_face_values(face_vals, num_fine_cells);
	malloc_star_values(star_vals, num_fine_cells);
	malloc_fluxes(fluxes, num_fine_cells);
	malloc_bar_values(bar_vals, num_fine_cells);
	malloc_scale_coefficients(scale_coeffs, totalScaleCoeffs);
	malloc_details(details, totalDetails);

	int* dry_cells = new int[num_fine_cells + 2];

	real* eta_temp = new real[num_fine_cells + 2];

	real* delta_west = new real[num_fine_cells + 1];
	real* delta_east = new real[num_fine_cells + 1];

	real dx_coarse = (sim_params.xmax - sim_params.xmin) / sim_params.cells;
	real dx_fine = dx_coarse / (1 << solver_params.L);

	real* dxFlattened = new real[totalScaleCoeffs];
	real* xFlattened = new real[totalScaleCoeffs];
	int* levelIndicesFlattened = new int[totalScaleCoeffs];

	real* norm_details = new real[totalDetails];
	int* sig_details = new int[totalDetails]();

	// allocate true/false buffer for dry cells
	bool* dry = new bool[num_fine_cells + 2]();

	// =========================================================== //

	get_nodal_values(nodal_vals, sim_params, num_fine_cells, bcs, dx_fine, test_case);

	real* xIntCoarse = new real[sim_params.cells + 1]();

	// coarse mesh
	for (i = 0; i < sim_params.cells + 1; i++)
	{
		xIntCoarse[i] = sim_params.xmin + dx_coarse * i;
	}

	// finest scale coefficients
	for (i = 0; i < num_fine_cells; i++)
	{
		assem_sol.q_BC[i + 1] = (nodal_vals.q[i] + nodal_vals.q[i + 1]) / 2;
		assem_sol.h_BC[i + 1] = (nodal_vals.h[i] + nodal_vals.h[i + 1]) / 2;
		assem_sol.z_BC[i + 1] = (nodal_vals.z[i] + nodal_vals.z[i + 1]) / 2;
	}

	// load the fine scale data
	for (int c = 0; c < sim_params.cells; c++)
	{
		int scaleStep = c * scaleCoeffsPerCell;

		for (int k = 0; k < (1 << solver_params.L); k++)
		{
			scale_coeffs.q[scaleStep + (1 << solver_params.L) - 1 + k] = assem_sol.q_BC[c * (1 << solver_params.L) + k + 1];
			scale_coeffs.eta[scaleStep + (1 << solver_params.L) - 1 + k] = assem_sol.h_BC[c * (1 << solver_params.L) + k + 1] + assem_sol.z_BC[c * (1 << solver_params.L) + k + 1];
			scale_coeffs.z[scaleStep + (1 << solver_params.L) - 1 + k] = assem_sol.z_BC[c * (1 << solver_params.L) + k + 1];
		}
	}

	while (timeNow < sim_params.simulationTime)
	{
		timeNow += dt;

		if (timeNow - sim_params.simulationTime > 0)
		{
			timeNow -= dt;
			dt = sim_params.simulationTime - timeNow;
			timeNow += dt;
		}

		// reset maxes
		qMax = 0;
		etaMax = 0;

		if (firstStep)
		{
			for (i = 0; i < num_fine_cells; i++)
			{
				qMax = max(qMax, abs(assem_sol.q_BC[i + 1]));
				zMax = max(zMax, abs(assem_sol.z_BC[i + 1]));
				etaMax = max(etaMax, abs(assem_sol.h_BC[i + 1] + assem_sol.z_BC[i + 1]));
			}
		}
		else
		{
			for (i = 0; i < assem_sol.length; i++)
			{
				real q = assem_sol.q_BC[i + 1];
				real eta = assem_sol.h_BC[i + 1] + assem_sol.z_BC[i + 1];

				scale_coeffs.q[assem_sol.activeIndices[i]] = q;
				scale_coeffs.eta[assem_sol.activeIndices[i]] = eta;

				qMax = max(qMax, abs(q));
				etaMax = max(etaMax, abs(eta));
			}
		}

		qMax = max(qMax, C(1.0));
		etaMax = max(etaMax, C(1.0));
		zMax = max(zMax, C(1.0));

		// thresholding details to zero for next step
		for (i = 0; i < totalDetails; i++)
		{
			details.q[i] = 0;
			details.eta[i] = 0;
		}


		// BEGIN ENCODING //

		for (int c = 0; c < sim_params.cells; c++)
		{
			int scaleStep = c * scaleCoeffsPerCell;
			int detailStep = c * detailsPerCell;

			for (int n = solver_params.L - 1; n >= 0; n--)
			{
				int currentLevStart = (1 << n) - 1;
				int currentLevEnd = (1 << (n + 1)) - 2;
				int kHigher = currentLevEnd + 1;

				for (int k = currentLevStart; k <= currentLevEnd; k++)
				{
					if (sig_details[detailStep + k] || firstStep)
					{
						real q1 = scale_coeffs.q[scaleStep + kHigher];
						real eta1 = scale_coeffs.eta[scaleStep + kHigher];

						real q2 = scale_coeffs.q[scaleStep + kHigher + 1];
						real eta2 = scale_coeffs.eta[scaleStep + kHigher + 1];

						scale_coeffs.q[scaleStep + k] = encodeScale(q1, q2);
						scale_coeffs.eta[scaleStep + k] = encodeScale(eta1, eta2);

						details.q[detailStep + k] = encodeDetail(q1, q2);
						details.eta[detailStep + k] = encodeDetail(eta1, eta2);

						if (firstStep)
						{
							real z1 = scale_coeffs.z[scaleStep + kHigher];
							real z2 = scale_coeffs.z[scaleStep + kHigher + 1];

							scale_coeffs.z[scaleStep + k] = encodeScale(z1, z2);
							details.z[detailStep + k] = encodeDetail(z1, z2);
						}
					}

					kHigher += 2;
				}
			}
		}

		// END ENCODING //

		firstStep = false;

		// zero significant details to reconstruct tree of details for the next iteration and normalise details
		for (i = 0; i < totalDetails; i++)
		{
			sig_details[i] = false;

			real a = abs(details.q[i]) / qMax;
			real b = abs(details.eta[i]) / etaMax;
			real c = abs(details.z[i]) / zMax;

			norm_details[i] = max(a, max(b, c));
		}

		// START PREDICTION //

		for (int c = 0; c < sim_params.cells; c++)
		{
			int scaleStep = c * scaleCoeffsPerCell;
			int detailStep = c * detailsPerCell;

			for (int n = 0; n < solver_params.L; n++)
			{
				int currentLevStart = (1 << n) - 1;
				int currentLevEnd = (1 << (n + 1)) - 2;

				for (int k = currentLevStart; k <= currentLevEnd; k++)
				{
					real epsilonLocal = solver_params.epsilon / ( 1 << (solver_params.L - n) );

					if (norm_details[detailStep + k] >= epsilonLocal)
					{
						sig_details[detailStep + k] = true;

						if (c + 1 < sim_params.cells && k == currentLevEnd) // if it's not the last cell and at the right-most edge
						{
							sig_details[(c + 1) * detailsPerCell + currentLevStart] = true; // make the subelement to the right cell also significant
						}

						if (c > 0 && k == currentLevStart) // if it's not the first cell and at the left-most edge
						{
							sig_details[(c - 1) * detailsPerCell + currentLevEnd] = true; // the subelement on the left cell is also significant
						}
					}
				}
			}
		}

		// END PREDICTION //

		

		// START REGULARISATION //

		for (int c = 0; c < sim_params.cells; c++)
		{
			int detailStep = c * detailsPerCell;

			for (int n = solver_params.L; n > 1; n--)
			{
				int k = (1 << (n - 1)) - 1;
				int kLower = (1 << (n - 2)) - 1;
				int currentLevEnd = (1 << n) - 2;

				for (k; k < currentLevEnd; k += 2)
				{
					if (sig_details[detailStep + k] || sig_details[detailStep + k + 1])
					{
						sig_details[detailStep + kLower] = true;
					}

					kLower++; // step along only the one parent element
				}
			}
		}

		// END REGULARISATION //

		

		// START EXTRA SIGNIFICANCE //

		for (int c = 0; c < sim_params.cells; c++)
		{
			int detailStep = c * detailsPerCell;

			for (int n = 0; n < solver_params.L; n++)
			{
				int currentLevStart = (1 << n) - 1;
				int currentLevEnd = (1 << (n + 1)) - 2;
				int kHigher = currentLevEnd + 1; // index of child element

				for (int k = currentLevStart; k <= currentLevEnd; k++)
				{
					if (sig_details[detailStep + k])
					{
						real mBar = 1.5;
						real epsilonLocal = solver_params.epsilon * pow(C(2.0), n - solver_params.L);

						if (norm_details[detailStep + k] >= epsilonLocal * pow(C(2.0), mBar + 1) && n + 1 != solver_params.L)
						{
							// if extra signficant child elements marked as active
							sig_details[detailStep + kHigher] = true;
							sig_details[detailStep + kHigher + 1] = true;
						}
					}

					kHigher += 2;
				}
			}
		}

		// END EXTRA SIGNIFICANCE //

		// reset since passing by reference
		assem_sol.length = 0;

		for (int c = 0; c < sim_params.cells; c++)
		{
			int scaleStep = c * scaleCoeffsPerCell;
			int detailStep = c * detailsPerCell;

			xFlattened[scaleStep] = (xIntCoarse[c] + xIntCoarse[c + 1]) / 2;
			dxFlattened[scaleStep] = dx_coarse;
			levelIndicesFlattened[scaleStep] = 0;

			// initially, n = k = 0
			treeTraversalDecode(solver_params, scale_coeffs, dxFlattened, xFlattened, levelIndicesFlattened, details,
				0, 0, detailStep, scaleStep, sig_details, assem_sol);
		}

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

		printf("Length: %d, step: %d, time step: %.15f, mass: %.1f, progress: %.17f%%\n", assem_sol.length, steps, dt, total_mass, timeNow / sim_params.simulationTime * 100);
	}

	// delete buffers

	free_nodal_values(nodal_vals);
	free_assembled_solution(assem_sol);
	free_face_values(face_vals);
	free_star_values(star_vals);
	free_bar_values(bar_vals);
	free_fluxes(fluxes);
	free_scale_coefficients(scale_coeffs);
	free_details(details);

	delete[] dry_cells;
	delete[] eta_temp;
	delete[] delta_west;
	delete[] delta_east;

	delete[] xIntCoarse;

	delete[] dxFlattened;
	delete[] xFlattened;
	delete[] levelIndicesFlattened;

	delete[] norm_details;
	delete[] sig_details;

	clock_t end = clock();

	real end_time = (real)(end - start) / CLOCKS_PER_SEC * C(1000.0);
	printf("Execution time measured using clock(): %f ms\n", end_time);

	return 0;
}