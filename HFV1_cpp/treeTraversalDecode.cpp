#include "treeTraversalDecode.h"

void treeTraversalDecode
(
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
)
{
	// initially, n = 0 and k = 0 i.e. the coarsest sub-element
	int currentLevStart = (1 << n) - 1;
	int currentLevEnd = (1 << (n + 1)) - 2;
	int nextLevStart = currentLevEnd + 1;

	// j is the local index along a level, whereas k is global along flattened tree
	int j = k - currentLevStart;

	int kHigher = 2 * j + nextLevStart;

	int level = levelIndicesFlattened[scaleStep + k];

	int detail = significantDetails[detailStep + k];

	if (level < solverParameters.L && detail)
	{
		// decode to generate multiresolution data
		real q1 = decode1(flattenedScaleCoeffs.q[scaleStep + k], flattenedDetails.q[detailStep + k]);
		real eta1 = decode1(flattenedScaleCoeffs.eta[scaleStep + k], flattenedDetails.eta[detailStep + k]);

		real q2 = decode2(flattenedScaleCoeffs.q[scaleStep + k], flattenedDetails.q[detailStep + k]);
		real eta2 = decode2(flattenedScaleCoeffs.eta[scaleStep + k], flattenedDetails.eta[detailStep + k]);

		// load data into the higher levels
		flattenedScaleCoeffs.q[scaleStep + kHigher] = q1;
		flattenedScaleCoeffs.eta[scaleStep + kHigher] = eta1;

		flattenedScaleCoeffs.q[scaleStep + kHigher + 1] = q2;
		flattenedScaleCoeffs.eta[scaleStep + kHigher + 1] = eta2;

		// generate multiresolution mesh coordinates
		real dxHigher = dxFlattened[scaleStep + k] / 2;
		dxFlattened[scaleStep + kHigher] = dxHigher;
		dxFlattened[scaleStep + kHigher + 1] = dxHigher;

		real x1 = xFlattened[scaleStep + k] - dxFlattened[scaleStep + kHigher] / 2;
		real x2 = xFlattened[scaleStep + k] + dxFlattened[scaleStep + kHigher] / 2;
		xFlattened[scaleStep + kHigher] = x1;
		xFlattened[scaleStep + kHigher + 1] = x2;

		// find index for sub-elements one level higher
		int nHigher = levelIndicesFlattened[scaleStep + k] + 1;
		levelIndicesFlattened[scaleStep + kHigher] = nHigher;
		levelIndicesFlattened[scaleStep + kHigher + 1] = nHigher;

		treeTraversalDecode
		(
			solverParameters, 
			flattenedScaleCoeffs, 
			dxFlattened, 
			xFlattened, 
			levelIndicesFlattened, 
			flattenedDetails, 
			n + 1, 
			kHigher, 
			detailStep, 
			scaleStep, 
			significantDetails, 
			assembledSolution
		);
		
		treeTraversalDecode
		(
			solverParameters,
			flattenedScaleCoeffs,
			dxFlattened,
			xFlattened,
			levelIndicesFlattened,
			flattenedDetails,
			n + 1,
			kHigher + 1,
			detailStep,
			scaleStep,
			significantDetails,
			assembledSolution
		);
	}
	else
	{
		// add the data to the array
		int idx = assembledSolution.length;

		real q = flattenedScaleCoeffs.q[scaleStep + k];
		real eta = flattenedScaleCoeffs.eta[scaleStep + k];
		real z = flattenedScaleCoeffs.z[scaleStep + k];
		real dxLocal = dxFlattened[scaleStep + k];
		real x = xFlattened[scaleStep + k];

		assembledSolution.q_BC[idx + 1] = q;
		assembledSolution.h_BC[idx + 1] = eta - z;
		assembledSolution.z_BC[idx + 1] = z;
		assembledSolution.dx_BC[idx + 1] = dxLocal;
		assembledSolution.x[idx + 1] = x;
		assembledSolution.activeIndices[idx] = scaleStep + k;
		assembledSolution.length++;
	}
}