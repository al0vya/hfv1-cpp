#pragma once

#include "real.h"
#include "SolverParameters.h"
#include "decode1.h"
#include "decode2.h"
#include "AssembledSolution.h"
#include "FlattenedDetails.h"
#include "FlattenedScaleCoeffs.h"

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