#include <stdlib.h>
#include <stdio.h>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <time.h>

#include "real.h"
#include "SimulationParameters.h"
#include "SolverParameters.h"
#include "BoundaryConditions.h"
#include "bedDataConservative.h"
#include "qInitial.h"
#include "hInitial.h"
#include "frictionImplicit.h"
#include "fluxHLL.h"
#include "encodeDetail.h"
#include "encodeScale.h"
#include "decode1.h"
#include "decode2.h"
#include "FlattenedScaleCoeffs.h"
#include "FlattenedDetails.h"
#include "AssembledSolution.h"
#include "treeTraversalDecode.h"

using namespace std;

void printTree(SolverParameters solverParameters, int* significantDetails, int scaleCoeffsPerCell, int detailsPerCell);
void printTreeOfReals(SolverParameters solverParameters, real* normalisedDetails, int scaleCoeffsPerCell, int detailsPerCell);

int main()
{
	clock_t start = clock();

	// quintessential for-loop index
	int i;
	int steps = 0;

	SimulationParameters simulationParameters;
	simulationParameters.xmin = 0;
	simulationParameters.xmax = 50;
	simulationParameters.simulationTime = C(100.0);
	simulationParameters.manning = C(0.02);

	SolverParameters solverParameters;
	solverParameters.cells = 1;
	solverParameters.CFL = C(0.33);
	solverParameters.tolDry = C(1e-3);
	solverParameters.g = C(9.80665);
	solverParameters.L = 10;

	BoundaryConditions bcs;
	bcs.hl = C(6.0);
	bcs.hr = C(0.0);

	bcs.ql = C(0.0);
	bcs.qr = C(0.0);

	bcs.reflectUp = C(1.0);
	bcs.reflectDown = C(1.0);

	bcs.hImposedUp = C(0.0);
	bcs.qxImposedUp = C(0.0);

	bcs.hImposedDown = C(0.0);
	bcs.qxImposedDown = C(0.0);

	real dxCoarse = (simulationParameters.xmax - simulationParameters.xmin) / solverParameters.cells;
	real dxFine = dxCoarse / (1 << solverParameters.L);

	solverParameters.epsilon = dxFine / 2;

	// number of cells/interfaces at finest resolution
	int cellsFine = solverParameters.cells * (1 << solverParameters.L);
	int interfacesFine = cellsFine + 1;

	real* xIntFine = new real[interfacesFine];
	real* qIntFine = new real[interfacesFine];
	real* hIntFine = new real[interfacesFine];
	real* zIntFine = new real[interfacesFine];

	// initialise finest mesh and flow nodes
	for (i = 0; i < interfacesFine; i++)
	{
		xIntFine[i] = simulationParameters.xmin + i * dxFine;
		zIntFine[i] = bedDataConservative(xIntFine[i]);
		hIntFine[i] = hInitial(bcs, zIntFine[i], xIntFine[i]);
		qIntFine[i] = qInitial(bcs, xIntFine[i]);
	}

	real* xIntCoarse = new real[solverParameters.cells + 1]();

	// coarse mesh
	for (i = 0; i < solverParameters.cells + 1; i++)
	{
		xIntCoarse[i] = simulationParameters.xmin + dxCoarse * i;
	}

	real* qScaleFine = new real[cellsFine];
	real* hScaleFine = new real[cellsFine];
	real* zScaleFine = new real[cellsFine];
	real* etaScaleFine = new real[cellsFine];

	// finest scale coefficients
	for (i = 0; i < cellsFine; i++)
	{
		qScaleFine[i] = (qIntFine[i] + qIntFine[i + 1]) / C(2.0);
		hScaleFine[i] = (hIntFine[i] + hIntFine[i + 1]) / 2;
		zScaleFine[i] = (zIntFine[i] + zIntFine[i + 1]) / 2;
		etaScaleFine[i] = zScaleFine[i] + hScaleFine[i];
	}

	int scaleCoeffsPerCell = (1 << (solverParameters.L + 1)) - 1;
	int totalScaleCoeffs = scaleCoeffsPerCell * solverParameters.cells;

	FlattenedScaleCoeffs flattenedScaleCoeffs;
	flattenedScaleCoeffs.q = new real[totalScaleCoeffs]();
	flattenedScaleCoeffs.eta = new real[totalScaleCoeffs]();
	flattenedScaleCoeffs.z = new real[totalScaleCoeffs]();

	int detailsPerCell = (1 << solverParameters.L) - 1;
	int totalDetails = detailsPerCell * solverParameters.cells;

	// load the fine scale data
	for (int c = 0; c < solverParameters.cells; c++)
	{
		int scaleStep = c * scaleCoeffsPerCell;

		for (int k = 0; k < (1 << solverParameters.L); k++)
		{
			flattenedScaleCoeffs.q[scaleStep + (1 << solverParameters.L) - 1 + k] = qScaleFine[c * (1 << solverParameters.L) + k];
			flattenedScaleCoeffs.eta[scaleStep + (1 << solverParameters.L) - 1 + k] = etaScaleFine[c * (1 << solverParameters.L) + k];
			flattenedScaleCoeffs.z[scaleStep + (1 << solverParameters.L) - 1 + k] = zScaleFine[c * (1 << solverParameters.L) + k];
		}
	}

	AssembledSolution assembledSolution;

	// defining to be the largest possible size for later reuse
	assembledSolution.qWithBC = new real[cellsFine + 2]();
	assembledSolution.hWithBC = new real[cellsFine + 2]();
	assembledSolution.zWithBC = new real[cellsFine + 2]();
	assembledSolution.dxLocalWithBC = new real[cellsFine + 2];
	assembledSolution.activeIndices = new int[cellsFine]();
	assembledSolution.length = 0;

	real* dxFlattened = new real[totalScaleCoeffs];
	real* xFlattened = new real[totalScaleCoeffs];
	int* levelIndicesFlattened = new int[totalScaleCoeffs];

	FlattenedDetails flattenedDetails;
	flattenedDetails.q = new real[totalDetails]();
	flattenedDetails.eta = new real[totalDetails]();
	flattenedDetails.z = new real[totalDetails]();

	real* normalisedDetails = new real[totalDetails];
	int* significantDetails = new int[totalDetails]();

	// allocate true/false buffer for dry cells
	bool* dry = new bool[cellsFine + 2]();

	real* etaTemp = new real[cellsFine + 2]();

	// allocating buffers for eastern and western interface values
	real* qEast = new real[cellsFine + 1];
	real* hEast = new real[cellsFine + 1];
	real* etaEast = new real[cellsFine + 1];

	real* qWest = new real[cellsFine + 1];
	real* hWest = new real[cellsFine + 1];
	real* etaWest = new real[cellsFine + 1];

	// allocating buffers for positivity preserving nodes
	real* qEastStar = new real[cellsFine + 1];
	real* hEastStar = new real[cellsFine + 1];

	real* qWestStar = new real[cellsFine + 1];
	real* hWestStar = new real[cellsFine + 1];

	real* zStarIntermediate = new real[cellsFine + 1];
	real* zStar = new real[cellsFine + 1];

	real* uWest = new real[cellsFine + 1];
	real* uEast = new real[cellsFine + 1];

	real* delta = new real[cellsFine + 1];

	// allocating buffers for numerical fluxes from HLL solver
	real* massFlux = new real[cellsFine + 1];
	real* momentumFlux = new real[cellsFine + 1];

	// allocating buffers for positivity preserving MODES
	real* hBar = new real[cellsFine];
	real* zBar = new real[cellsFine];

	bool firstStep = true;
	real timeNow = 0;
	real dt = C(1e-4);

	while (timeNow < simulationParameters.simulationTime)
	{
		timeNow += dt;

		if (timeNow - simulationParameters.simulationTime > 0)
		{
			timeNow -= dt;
			dt = simulationParameters.simulationTime - timeNow;
			timeNow += dt;
		}

		// reset maxes
		real qMax = 0;
		real zMax = 0;
		real etaMax = 0;

		if (firstStep)
		{
			for (i = 0; i < cellsFine; i++)
			{
				qMax = max(qMax, abs(qScaleFine[i]));
				zMax = max(zMax, abs(zScaleFine[i]));
				etaMax = max(etaMax, abs(etaScaleFine[i]));
			}
		}
		else
		{
			for (i = 0; i < assembledSolution.length; i++)
			{
				real q = assembledSolution.qWithBC[i + 1];
				real eta = assembledSolution.hWithBC[i + 1] + assembledSolution.zWithBC[i + 1];
				real z = assembledSolution.zWithBC[i + 1];

				flattenedScaleCoeffs.q[assembledSolution.activeIndices[i]] = q;
				flattenedScaleCoeffs.eta[assembledSolution.activeIndices[i]] = eta;
				flattenedScaleCoeffs.z[assembledSolution.activeIndices[i]] = z;

				qMax = max(qMax, abs(q));
				etaMax = max(etaMax, abs(eta));
				zMax = max(zMax, abs(z));
			}
		}

		// thresholding details to zero for next step
		for (i = 0; i < totalDetails; i++)
		{
			flattenedDetails.q[i] = 0;
			flattenedDetails.eta[i] = 0;
			flattenedDetails.z[i] = 0;
		}


		// BEGIN ENCODING //

		for (int c = 0; c < solverParameters.cells; c++)
		{
			int scaleStep = c * scaleCoeffsPerCell;
			int detailStep = c * detailsPerCell;

			for (int n = solverParameters.L - 1; n >= 0; n--)
			{
				int currentLevStart = (1 << n) - 1;
				int currentLevEnd = (1 << (n + 1)) - 2;
				int kHigher = currentLevEnd + 1;

				for (int k = currentLevStart; k <= currentLevEnd; k++)
				{
					if (significantDetails[detailStep + k] || firstStep)
					{
						real q1 = flattenedScaleCoeffs.q[scaleStep + kHigher];
						real eta1 = flattenedScaleCoeffs.eta[scaleStep + kHigher];
						real z1 = flattenedScaleCoeffs.z[scaleStep + kHigher];

						real q2 = flattenedScaleCoeffs.q[scaleStep + kHigher + 1];
						real eta2 = flattenedScaleCoeffs.eta[scaleStep + kHigher + 1];
						real z2 = flattenedScaleCoeffs.z[scaleStep + kHigher + 1];

						flattenedScaleCoeffs.q[scaleStep + k] = encodeScale(q1, q2);
						flattenedScaleCoeffs.eta[scaleStep + k] = encodeScale(eta1, eta2);
						flattenedScaleCoeffs.z[scaleStep + k] = encodeScale(z1, z2);

						flattenedDetails.q[detailStep + k] = encodeDetail(q1, q2);
						flattenedDetails.eta[detailStep + k] = encodeDetail(eta1, eta2);
						flattenedDetails.z[detailStep + k] = encodeDetail(z1, z2);
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
			significantDetails[i] = false;

			real a = (qMax != 0) ? flattenedDetails.q[i] / qMax : 0;
			real b = (etaMax != 0) ? flattenedDetails.eta[i] / etaMax : 0;
			real c = (zMax != 0) ? flattenedDetails.z[i] / zMax : 0;

			normalisedDetails[i] = max(a, max(b, c));
		}

		//printTreeOfReals(solverParameters, normalisedDetails, scaleCoeffsPerCell, detailsPerCell);

		// START PREDICTION //

		for (int c = 0; c < solverParameters.cells; c++)
		{
			int scaleStep = c * scaleCoeffsPerCell;
			int detailStep = c * detailsPerCell;

			for (int n = 0; n < solverParameters.L; n++)
			{
				int currentLevStart = (1 << n) - 1;
				int currentLevEnd = (1 << (n + 1)) - 2;

				for (int k = currentLevStart; k <= currentLevEnd; k++)
				{
					real epsilonLocal = solverParameters.epsilon * pow(C(2.0), n - solverParameters.L);

					if (normalisedDetails[detailStep + k] > epsilonLocal)
					{
						significantDetails[detailStep + k] = true;

						if (c + 1 < solverParameters.cells && k == currentLevEnd) // if it's not the last cell and at the right-most edge
						{
							significantDetails[(c + 1) * detailsPerCell + currentLevStart] = true; // make the subelement to the right cell also significant
						}

						if (c > 0 && k == currentLevStart) // if it's not the first cell and at the left-most edge
						{
							significantDetails[(c - 1) * detailsPerCell + currentLevEnd] = true; // the subelement on the left cell is also significant
						}
					}
				}
			}
		}

		// END PREDICTION //

		printTree(solverParameters, significantDetails, scaleCoeffsPerCell, detailsPerCell);

		// START REGULARISATION //

		for (int c = 0; c < solverParameters.cells; c++)
		{
			int detailStep = c * detailsPerCell;

			for (int n = solverParameters.L; n > 1; n--)
			{
				int k = (1 << (n - 1)) - 1;
				int kLower = (1 << (n - 2)) - 1;
				int currentLevEnd = (1 << n) - 2;

				for (k; k < currentLevEnd; k += 2)
				{
					if (significantDetails[detailStep + k] || significantDetails[detailStep + k + 1])
					{
						significantDetails[detailStep + kLower] = true;
					}

					kLower++; // step along only the one parent element
				}
			}
		}

		// END REGULARISATION //

		// START EXTRA SIGNIFICANCE //

		for (int c = 0; c < solverParameters.cells; c++)
		{
			int detailStep = c * detailsPerCell;

			for (int n = 0; n < solverParameters.L; n++)
			{
				int currentLevStart = (1 << n) - 1;
				int currentLevEnd = (1 << (n + 1)) - 2;
				int kHigher = currentLevEnd + 1; // index of child element

				for (int k = currentLevStart; k <= currentLevEnd; k++)
				{
					if (significantDetails[detailStep + k])
					{
						real mBar = 1.5;
						real epsilonLocal = solverParameters.epsilon * pow(C(2.0), n - solverParameters.L);

						if (normalisedDetails[detailStep + k] >= epsilonLocal * pow(2, mBar + 1) && n + 1 != solverParameters.L)
						{
							// if extra signficant child elements marked as active
							significantDetails[detailStep + kHigher] = true;
							significantDetails[detailStep + kHigher + 1] = true;
						}
					}

					kHigher += 2;
				}
			}
		}

		// END EXTRA SIGNIFICANCE //

		// reset since passing by reference
		assembledSolution.length = 0;

		for (int c = 0; c < solverParameters.cells; c++)
		{
			int scaleStep = c * scaleCoeffsPerCell;
			int detailStep = c * detailsPerCell;

			xFlattened[scaleStep] = (xIntCoarse[c] + xIntCoarse[c + 1]) / 2;
			dxFlattened[scaleStep] = dxCoarse;
			levelIndicesFlattened[scaleStep] = 0;

			// initially, n = k = 0
			treeTraversalDecode(solverParameters, flattenedScaleCoeffs, dxFlattened, xFlattened, levelIndicesFlattened, flattenedDetails,
				0, 0, detailStep, scaleStep, significantDetails, assembledSolution);
		}

		// adding ghost boundary conditions
		int end = assembledSolution.length + 1;

		assembledSolution.qWithBC[0] = (bcs.qxImposedUp > 0) ? bcs.qxImposedUp : bcs.reflectUp * assembledSolution.qWithBC[1];
		assembledSolution.qWithBC[end] = (bcs.qxImposedDown > 0) ? bcs.qxImposedDown : bcs.reflectDown * assembledSolution.qWithBC[end - 1]; // recall there are cells + 2 elements inc BCs

		assembledSolution.hWithBC[0] = (bcs.hImposedDown > 0) ? bcs.hImposedUp : assembledSolution.hWithBC[1];
		assembledSolution.hWithBC[end] = (bcs.hImposedDown > 0) ? bcs.hImposedDown : assembledSolution.hWithBC[end - 1];

		assembledSolution.zWithBC[0] = assembledSolution.zWithBC[1];
		assembledSolution.zWithBC[end] = assembledSolution.zWithBC[end - 1];

		assembledSolution.dxLocalWithBC[0] = assembledSolution.dxLocalWithBC[1];
		assembledSolution.dxLocalWithBC[end] = assembledSolution.dxLocalWithBC[end - 1];

		// extract upwind and downwind modes
		real hWestUpwind = assembledSolution.hWithBC[0];
		real hEastDownwind = assembledSolution.hWithBC[end];

		if (simulationParameters.manning > 0)
		{
			for (i = 1; i < assembledSolution.length + 1; i++)
			{
				assembledSolution.qWithBC[i] += frictionImplicit(simulationParameters, solverParameters, dt, assembledSolution.hWithBC[i], assembledSolution.qWithBC[i]);
			}
		}

		// initialising dry vs wet cells and etaTemp, ignore ghost cells
		for (i = 1; i < assembledSolution.length + 1; i++)
		{
			real hLocal = assembledSolution.hWithBC[i];
			real hBackward = assembledSolution.hWithBC[i - 1];
			real hForward = assembledSolution.hWithBC[i + 1];

			real hMax = max(hLocal, hBackward);
			hMax = max(hForward, hMax);

			dry[i] = (hMax <= solverParameters.tolDry);

			etaTemp[i] = assembledSolution.hWithBC[i] + assembledSolution.zWithBC[i];
			int dummy = 1;
		}

		// initialising interface values
		for (i = 0; i < assembledSolution.length + 1; i++)
		{
			qEast[i] = assembledSolution.qWithBC[i + 1];
			hEast[i] = assembledSolution.hWithBC[i + 1];
			etaEast[i] = etaTemp[i + 1];

			qWest[i] = assembledSolution.qWithBC[i];
			hWest[i] = assembledSolution.hWithBC[i];
			etaWest[i] = etaTemp[i];
		}

		// correcting downwind and upwind eta values
		etaEast[assembledSolution.length] = etaTemp[assembledSolution.length] - assembledSolution.hWithBC[assembledSolution.length] + hEastDownwind;
		etaWest[0] = etaTemp[1] - assembledSolution.hWithBC[1] + hWestUpwind;

		for (int i = 0; i < assembledSolution.length + 1; i++)
		{
			// initialising velocity interface values
			uWest[i] = (hWest[i] <= solverParameters.tolDry) ? 0 : qWest[i] / hWest[i];
			uEast[i] = (hEast[i] <= solverParameters.tolDry) ? 0 : qEast[i] / hEast[i];

			// intermediate calculations
			real a = etaWest[i] - hWest[i];
			real b = etaEast[i] - hEast[i];
			zStarIntermediate[i] = max(a, b);
			a = etaWest[i] - zStarIntermediate[i];
			b = etaEast[i] - zStarIntermediate[i];

			// positivity-preserving nodes
			hWestStar[i] = max(C(0.0), a);
			hEastStar[i] = max(C(0.0), b);

			delta[i] = max(C(0.0), -a) + max(C(0.0), -b);

			qWestStar[i] = uWest[i] * hWestStar[i];
			qEastStar[i] = uEast[i] * hEastStar[i];

			zStar[i] = zStarIntermediate[i] - delta[i];
		}

		// initialising numerical fluxes
		fluxHLL(assembledSolution.length, solverParameters, hWestStar, hEastStar, qWestStar, qEastStar, uWest, uEast, massFlux, momentumFlux);

		for (int i = 0; i < assembledSolution.length; i++)
		{
			// essentially 0th order projection but taking into account east/west locality
			hBar[i] = (hEastStar[i] + hWestStar[i + 1]) / 2;

			// 1st order projection
			zBar[i] = (zStar[i + 1] - zStar[i]) / (2 * sqrt(C(3.0)));
		}

		// FV1 operator increment, skip ghosts cells
		for (i = 1; i < assembledSolution.length + 1; i++)
		{
			// skip increment in dry cells
			if (!dry[i])
			{
				real massIncrement = -(1 / assembledSolution.dxLocalWithBC[i]) * (massFlux[i] - massFlux[i - 1]);
				real momentumIncrement = -(1 / assembledSolution.dxLocalWithBC[i]) * (momentumFlux[i] - momentumFlux[i - 1] + 2 * sqrt(C(3.0)) * solverParameters.g * hBar[i - 1] * zBar[i - 1]);

				assembledSolution.hWithBC[i] += dt * massIncrement;
				assembledSolution.qWithBC[i] = (assembledSolution.hWithBC[i] > solverParameters.tolDry) ? assembledSolution.qWithBC[i] + dt * momentumIncrement : 0;
			}
		}

		// CFL time step adjustment
		dt = 1e9;
		real totalMass = 0;
		steps++;

		for (i = 1; i < assembledSolution.length + 1; i++)
		{
			if (assembledSolution.hWithBC[i] > solverParameters.tolDry)
			{
				real u = assembledSolution.qWithBC[i] / assembledSolution.hWithBC[i];
				real dtCFL = solverParameters.CFL * assembledSolution.dxLocalWithBC[i] / (abs(u) + sqrt(solverParameters.g * assembledSolution.hWithBC[i]));
				dt = min(dt, dtCFL);
			}

			totalMass += assembledSolution.hWithBC[i] * assembledSolution.dxLocalWithBC[i];
		}

	/*	for (i = 1; i < assembledSolution.length + 1; i++)
		{
			printf("%0.2f, ", assembledSolution.hWithBC[i]);
		}
		printf("\n");*/
		//printf("Mass: %.17g, dt: %f, simulation time: %f s\n", totalMass, dt, timeNow);
		printf("%d\n", assembledSolution.length);
	}

	// delete buffers
	delete[] xIntFine;
	delete[] qIntFine;
	delete[] hIntFine;
	delete[] zIntFine;

	delete[] qScaleFine;
	delete[] hScaleFine;
	delete[] zScaleFine;
	delete[] etaScaleFine;

	delete[] xIntCoarse;

	delete[] flattenedScaleCoeffs.q;
	delete[] flattenedScaleCoeffs.eta;
	delete[] flattenedScaleCoeffs.z;

	delete[] flattenedDetails.q;
	delete[] flattenedDetails.eta;
	delete[] flattenedDetails.z;

	delete[] assembledSolution.qWithBC;
	delete[] assembledSolution.hWithBC;
	delete[] assembledSolution.zWithBC;
	delete[] assembledSolution.dxLocalWithBC;
	delete[] assembledSolution.activeIndices;

	delete[] dxFlattened;
	delete[] xFlattened;
	delete[] levelIndicesFlattened;

	delete[] normalisedDetails;
	delete[] significantDetails;

	delete[] dry;

	delete[] etaTemp;

	delete[] qEast;
	delete[] hEast;
	delete[] etaEast;

	delete[] qWest;
	delete[] hWest;
	delete[] etaWest;

	delete[] qWestStar;
	delete[] hWestStar;

	delete[] qEastStar;
	delete[] hEastStar;

	delete[] zStarIntermediate;
	delete[] zStar;

	delete[] uWest;
	delete[] uEast;

	delete[] delta;

	delete[] massFlux;
	delete[] momentumFlux;

	delete[] hBar;
	delete[] zBar;

	clock_t end = clock();

	real time = (real)(end - start) / CLOCKS_PER_SEC * C(1000.0);
	printf("Execution time measured using clock(): %f ms\n", time);

	return 0;
}

void printTree(SolverParameters solverParameters, int* significantDetails, int scaleCoeffsPerCell, int detailsPerCell)
{
	for (int c = 0; c < solverParameters.cells; c++)
	{
		int scaleStep = c * scaleCoeffsPerCell;
		int detailStep = c * detailsPerCell;

		for (int n = solverParameters.L - 1; n >= 0; n--)
		{
			int currentLevStart = (1 << n) - 1;
			int currentLevEnd = (1 << (n + 1)) - 2;

			for (int k = currentLevStart; k <= currentLevEnd; k++)
			{
				printf("| %d |", significantDetails[detailStep + k]);
			}

			printf("\n");
		}

		printf("\n");
	}
}

void printTreeOfReals(SolverParameters solverParameters, real* normalisedDetails, int scaleCoeffsPerCell, int detailsPerCell)
{
	for (int c = 0; c < solverParameters.cells; c++)
	{
		int detailStep = c * detailsPerCell;

		for (int n = solverParameters.L - 1; n >= 0; n--)
		{
			int currentLevStart = (1 << n) - 1;
			int currentLevEnd = (1 << (n + 1)) - 2;

			for (int k = currentLevStart; k <= currentLevEnd; k++)
			{
				printf("| %0.2f |", normalisedDetails[detailStep + k]);
			}

			printf("\n");
		}

		printf("\n");
	}
}