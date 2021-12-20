#pragma once

#include "../classes/real.h"
#include "../classes/SolverParameters.h"
#include "decode1.h"
#include "decode2.h"
#include "../classes/AssembledSolution.h"
#include "../classes/FlattenedDetails.h"
#include "../classes/FlattenedScaleCoeffs.h"

void traverse_tree_decode(
	SolverParameters solverParameters,

	FlattenedScaleCoeffs flattenedScaleCoeffs,

	real* dxFlattened,
	real* xFlattened,
	int* levelIndicesFlattened,

	FlattenedDetails flattenedDetails,

	int n,
	int k,

	int detailStep,
	int scaleStep,

	int* significantDetails,

	AssembledSolution& assembledSolution
);