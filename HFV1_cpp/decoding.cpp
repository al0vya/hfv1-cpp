#include "decoding.h"

void decoding
(
	AssembledSolution&    assem_sol, 
	SimulationParameters& sim_params, 
	int&                  scale_coeffs_per_cell, 
	int&                  details_per_cell,
	real*&                xFlattened, 
	real*&                xIntCoarse, 
	real*&                dxFlattened, 
	real&                 dx_coarse, 
	int*&                 levelIndicesFlattened, 
	int*&                 sig_details, 
	SolverParameters&     solver_params, 
	FlattenedScaleCoeffs  scale_coeffs, 
	FlattenedDetails      details
)
{
	// reset since passing by reference
	assem_sol.length = 0;

	for (int cell = 0; cell < sim_params.cells; cell++)
	{
		int scaleStep = cell * scale_coeffs_per_cell;
		int detailStep = cell * details_per_cell;

		xFlattened[scaleStep] = (xIntCoarse[cell] + xIntCoarse[cell + 1]) / 2;
		dxFlattened[scaleStep] = dx_coarse;
		levelIndicesFlattened[scaleStep] = 0;

		// initially, n = k = 0
		treeTraversalDecode
		(
			solver_params,
			scale_coeffs,
			dxFlattened,
			xFlattened,
			levelIndicesFlattened,
			details,
			0,
			0,
			detailStep,
			scaleStep,
			sig_details,
			assem_sol
		);
	}
}