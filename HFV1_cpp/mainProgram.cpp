#include <stdlib.h>
#include <stdio.h>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <fstream>
#include <time.h>

#include "real.h"
#include "SimulationParameters.h"
#include "SolverParameters.h"
#include "BoundaryConditions.h"
#include "set_simulation_parameters_for_test_case.h"
#include "set_solver_parameters_for_test_case.h"
#include "set_boundary_conditions_for_test_case.h"
#include "bedDataConservative.h"
#include "hInitialConservative.h"
#include "hInitialOvertopping.h"
#include "qInitial.h"
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
#include "writeSolutionToFile.h"

using namespace std;

void printTree(SimulationParameters simulationParameters, SolverParameters solverParameters, int* significantDetails, int detailsPerCell);
void printTreeOfReals(SimulationParameters simulationParameters, SolverParameters solverParameters, real* details, int detailsPerCell);
void printScaleCoeffs(SimulationParameters simulationParameters, SolverParameters solverParameters, real* scaleCoeffs, int scalesCoeffsPerCell);

int main()
{
	clock_t start = clock();

	// quintessential for-loop index
	int i;
	int steps = 0;

	std::cout << "Please enter a number between 1 and 6 to select a test case.\n"
		"1: Wet dam break\n"
		"2: Dry dam break\n"
		"3: Dry dam break with friction\n"
		"4: Wet lake-at-rest (C-property)\n"
		"5: Wet/dry lake-at-rest\n"
		"6: Building overtopping\n";

	int test_case_selection;

	std::cin >> test_case_selection;

	if (!std::cin || test_case_selection > 6 || test_case_selection < 1)
	{
		std::cout << "Error: please rerun and enter a number between 1 and 6. Exiting program.\n";

		return -1;
	}

	std::cout << "Please enter the number of mother elements.\n";

	int number_of_cells;

	std::cin >> number_of_cells;

	if (!std::cin || number_of_cells < 1)
	{
		std::cout << "Error: please rerun and enter a integer value. Exiting program.\n";

		return -1;
	}

	std::cout << "Please select a maximum refinement level.\n";

	int user_input_max_refinement_level;

	std::cin >> user_input_max_refinement_level;

	if (!std::cin || user_input_max_refinement_level < 1)
	{
		std::cout << "Error: please rerun and enter a integer value. Exiting program.\n";

		return -1;
	}

	std::cout << "Please enter an error threshold.\n";

	real user_input_epsilon;

	std::cin >> user_input_epsilon;

	if (!std::cin)
	{
		std::cout << "Error: please rerun and enter a float or double value. Exiting program.\n";

		return -1;
	}

	SimulationParameters simulationParameters = set_simulation_parameters_for_test_case(test_case_selection, number_of_cells);
	SolverParameters solverParameters = set_solver_parameters_for_test_case(user_input_max_refinement_level, user_input_epsilon);
	BoundaryConditions bcs = set_boundary_conditions_for_test_case(test_case_selection);

	real dxCoarse = (simulationParameters.xmax - simulationParameters.xmin) / simulationParameters.cells;
	real dxFine = dxCoarse / (1 << solverParameters.L);

	// number of cells/interfaces at finest resolution
	int cellsFine = simulationParameters.cells * (1 << solverParameters.L);
	int interfacesFine = cellsFine + 1;

	real* xIntFine = new real[interfacesFine];
	real* qIntFine = new real[interfacesFine];
	real* hIntFine = new real[interfacesFine];
	real* zIntFine = new real[interfacesFine];

	// initialise finest mesh and flow nodes
	for (i = 0; i < interfacesFine; i++)
	{
		xIntFine[i] = simulationParameters.xmin + i * dxFine;
	}

	switch (test_case_selection)
	{
	case 1:
	case 2:
	case 3:
		for (i = 0; i < interfacesFine; i++)
		{
			zIntFine[i] = 0;
			hIntFine[i] = hInitialOvertopping(bcs, zIntFine[i], xIntFine[i]);
			qIntFine[i] = qInitial(bcs, xIntFine[i]);
		}
		break;
	case 4:
	case 5:
		for (i = 0; i < interfacesFine; i++)
		{
			zIntFine[i] = bedDataConservative(xIntFine[i]);
			hIntFine[i] = hInitialConservative(bcs, zIntFine[i], xIntFine[i]);
			qIntFine[i] = qInitial(bcs, xIntFine[i]);
		}
		break;
	case 6:
		for (i = 0; i < interfacesFine; i++)
		{
			zIntFine[i] = bedDataConservative(xIntFine[i]);
			hIntFine[i] = hInitialOvertopping(bcs, zIntFine[i], xIntFine[i]);
			qIntFine[i] = qInitial(bcs, xIntFine[i]);
		}
		break;
	default:
		break;
	}

	real* xIntCoarse = new real[simulationParameters.cells + 1]();

	// coarse mesh
	for (i = 0; i < simulationParameters.cells + 1; i++)
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
		qScaleFine[i] = (qIntFine[i] + qIntFine[i + 1]) / 2;
		hScaleFine[i] = (hIntFine[i] + hIntFine[i + 1]) / 2;
		zScaleFine[i] = (zIntFine[i] + zIntFine[i + 1]) / 2;
		etaScaleFine[i] = zScaleFine[i] + hScaleFine[i];
	}

	int scaleCoeffsPerCell = (1 << (solverParameters.L + 1)) - 1;
	int totalScaleCoeffs = scaleCoeffsPerCell * simulationParameters.cells;

	FlattenedScaleCoeffs flattenedScaleCoeffs;
	flattenedScaleCoeffs.q = new real[totalScaleCoeffs]();
	flattenedScaleCoeffs.eta = new real[totalScaleCoeffs]();
	flattenedScaleCoeffs.z = new real[totalScaleCoeffs]();

	int detailsPerCell = (1 << solverParameters.L) - 1;
	int totalDetails = detailsPerCell * simulationParameters.cells;

	// load the fine scale data
	for (int c = 0; c < simulationParameters.cells; c++)
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
	assembledSolution.x = new real[cellsFine + 2];
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

	real* zWestStar = new real[cellsFine + 1];
	real* zEastStar = new real[cellsFine + 1];

	real* deltaWest = new real[cellsFine + 1];
	real* deltaEast = new real[cellsFine + 1];

	// allocating buffers for numerical fluxes from HLL solver
	real* massFlux = new real[cellsFine + 1];
	real* momentumFlux = new real[cellsFine + 1];

	// allocating buffers for positivity preserving MODES
	real* hBar = new real[cellsFine];
	real* zBar = new real[cellsFine];

	real qMax = 0;
	real zMax = 0;
	real etaMax = 0;

	bool firstStep = true;
	real timeNow = 0;
	real dt = C(1e-4);

	ofstream data;

	data.open("clock_time_vs_sim_time.csv");

	data.precision(15);

	data << "sim_time,clock_time" << std::endl;

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
		qMax = 0;
		etaMax = 0;

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

				flattenedScaleCoeffs.q[assembledSolution.activeIndices[i]] = q;
				flattenedScaleCoeffs.eta[assembledSolution.activeIndices[i]] = eta;

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
			flattenedDetails.q[i] = 0;
			flattenedDetails.eta[i] = 0;
		}


		// BEGIN ENCODING //

		for (int c = 0; c < simulationParameters.cells; c++)
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

						real q2 = flattenedScaleCoeffs.q[scaleStep + kHigher + 1];
						real eta2 = flattenedScaleCoeffs.eta[scaleStep + kHigher + 1];

						flattenedScaleCoeffs.q[scaleStep + k] = encodeScale(q1, q2);
						flattenedScaleCoeffs.eta[scaleStep + k] = encodeScale(eta1, eta2);

						flattenedDetails.q[detailStep + k] = encodeDetail(q1, q2);
						flattenedDetails.eta[detailStep + k] = encodeDetail(eta1, eta2);

						if (firstStep)
						{
							real z1 = flattenedScaleCoeffs.z[scaleStep + kHigher];
							real z2 = flattenedScaleCoeffs.z[scaleStep + kHigher + 1];

							flattenedScaleCoeffs.z[scaleStep + k] = encodeScale(z1, z2);
							flattenedDetails.z[detailStep + k] = encodeDetail(z1, z2);
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
			significantDetails[i] = false;

			real a = abs(flattenedDetails.q[i]) / qMax;
			real b = abs(flattenedDetails.eta[i]) / etaMax;
			real c = abs(flattenedDetails.z[i]) / zMax;

			normalisedDetails[i] = max(a, max(b, c));
		}

		// START PREDICTION //

		for (int c = 0; c < simulationParameters.cells; c++)
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

						if (c + 1 < simulationParameters.cells && k == currentLevEnd) // if it's not the last cell and at the right-most edge
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


		// START REGULARISATION //

		for (int c = 0; c < simulationParameters.cells; c++)
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

		for (int c = 0; c < simulationParameters.cells; c++)
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

						if (normalisedDetails[detailStep + k] >= epsilonLocal * pow(C(2.0), mBar + 1) && n + 1 != solverParameters.L)
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

		for (int c = 0; c < simulationParameters.cells; c++)
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
			real uWest = (hWest[i] <= solverParameters.tolDry) ? 0 : qWest[i] / hWest[i];
			real uEast = (hEast[i] <= solverParameters.tolDry) ? 0 : qEast[i] / hEast[i];

			// intermediate calculations
			real zWest = etaWest[i] - hWest[i];
			real zEast = etaEast[i] - hEast[i];

			real zIntermediate = max(zWest, zEast);

			deltaWest[i] = max(C(0.0), -(etaWest[i] - zIntermediate));
			deltaEast[i] = max(C(0.0), -(etaEast[i] - zIntermediate));

			// positivity-preserving nodes
			hWestStar[i] = max(C(0.0), etaWest[i] - zIntermediate);
			qWestStar[i] = uWest * hWestStar[i];

			hEastStar[i] = max(C(0.0), etaEast[i] - zIntermediate);
			qEastStar[i] = uEast * hEastStar[i];

			zWestStar[i] = zIntermediate - deltaWest[i];
			zEastStar[i] = zIntermediate - deltaEast[i];
		}

		// initialising numerical fluxes
		fluxHLL(assembledSolution, solverParameters, hWestStar, hEastStar, qWestStar, qEastStar, massFlux, momentumFlux);

		for (int i = 0; i < assembledSolution.length; i++)
		{
			// essentially 0th order projection but taking into account east/west locality
			hBar[i] = (hWestStar[i + 1] + hEastStar[i]) / 2;

			// 1st order projection
			zBar[i] = (zWestStar[i + 1] - zEastStar[i]) / (2 * sqrt(C(3.0)));
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
				assembledSolution.qWithBC[i] = (assembledSolution.hWithBC[i] <= solverParameters.tolDry) ? 0 : assembledSolution.qWithBC[i] + dt * momentumIncrement;
			}
		}

		// CFL time step adjustment
		dt = 1e9;
		real totalMass = 0;

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

		steps++;

		if (steps % 10 == 0)
		{
			clock_t current = clock();
			real time = (real)(current - start) / CLOCKS_PER_SEC * C(1000.0);

			data << timeNow << "," << time << endl;
		}

		printf("Length: %d, step: %d, time step: %.15f, mass: %.15f, progress: %.17f%%\n", assembledSolution.length, steps, dt, totalMass, timeNow / simulationParameters.simulationTime * 100);
	}

	clock_t current = clock();
	real time = (real)(current - start) / CLOCKS_PER_SEC * C(1000.0);

	data << timeNow << "," << time << endl;

	data.close();

	writeSolutionToFile(assembledSolution);

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
	delete[] assembledSolution.x;
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

	delete[] zWestStar;
	delete[] zEastStar;

	delete[] deltaWest;
	delete[] deltaEast;

	delete[] massFlux;
	delete[] momentumFlux;

	delete[] hBar;
	delete[] zBar;

	clock_t end = clock();

	real end_time = (real)(end - start) / CLOCKS_PER_SEC * C(1000.0);
	printf("Execution time measured using clock(): %f ms\n", end_time);

	return 0;
}

void printTree(SimulationParameters simulationParameters, SolverParameters solverParameters, int* significantDetails, int detailsPerCell)
{
	for (int c = 0; c < simulationParameters.cells; c++)
	{
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

void printTreeOfReals(SimulationParameters simulationParameters, SolverParameters solverParameters, real* details, int detailsPerCell)
{
	for (int c = 0; c < simulationParameters.cells; c++)
	{
		int detailStep = c * detailsPerCell;

		for (int n = solverParameters.L - 1; n >= 0; n--)
		{
			int currentLevStart = (1 << n) - 1;
			int currentLevEnd = (1 << (n + 1)) - 2;

			for (int k = currentLevStart; k <= currentLevEnd; k++)
			{
				printf("| %.17f, %d |", (details[detailStep + k]), k);
			}

			printf("\n");
		}

		printf("\n");
	}
}

void printScaleCoeffs(SimulationParameters simulationParameters, SolverParameters solverParameters, real* scaleCoeffs, int scaleCoeffsPerCell)
{
	for (int c = 0; c < simulationParameters.cells; c++)
	{
		int scaleStep = c * scaleCoeffsPerCell;

		for (int n = solverParameters.L; n >= 0; n--)
		{
			int currentLevStart = (1 << n) - 1;
			int currentLevEnd = (1 << (n + 1)) - 2;

			for (int k = currentLevStart; k <= currentLevEnd; k++)
			{
				printf("| %.20g |", scaleCoeffs[scaleStep + k]);
			}

			printf("\n");
		}

		printf("\n");
	}
}