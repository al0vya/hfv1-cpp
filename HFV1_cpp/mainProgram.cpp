#include <stdlib.h>
#include <stdio.h>
#include <algorithm>
#include <cmath>
#include <vector>
#include <iostream>

#include "real.h"
#include "SimulationParameters.h"
#include "SolverParameters.h"
#include "BoundaryConditions.h"
#include "bedDataConservative.h"
#include "qInitial.h"
#include "hInitial.h"
#include "frictionImplicit.h"
#include "uFaceValues.h"
#include "fluxHLL.h"
#include "myPowInt.h"
#include "encodeDetail.h"
#include "encodeScale.h"
#include "decode1.h"
#include "decode2.h"

using namespace std;

template <typename A, typename B>
void zip(const vector<A>& a, const vector<B>& b, vector<pair<A, B>>& zipped);
template <typename A, typename B>
void unzip(const vector<pair<A, B>>& zipped, vector<A>& a, vector<B>& b);

int main()
{
	// quintessential for-loop index
	int i;

	SimulationParameters simulationParameters;
	simulationParameters.cells = 10;
	simulationParameters.xmin = 0;
	simulationParameters.xmax = 50;
	simulationParameters.simulationTime = C(500.0);
	simulationParameters.manning = C(0.02);
	
	SolverParameters solverParameters;
	solverParameters.CFL = C(0.33);
	solverParameters.tolDry = C(1e-3);
	solverParameters.g = C(9.80665);
	solverParameters.L = 5;

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

	real dxCoarse = (simulationParameters.xmax - simulationParameters.xmin) / simulationParameters.cells;
	real dxFine = dxCoarse / myPowInt(2, solverParameters.L);

	solverParameters.epsilon = dxFine / 2;

	// number of cells/interfaces at finest resolution
	int cellsFine = simulationParameters.cells * myPowInt(2, solverParameters.L);
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

	real* qScaleFine = new real[cellsFine];
	real* hScaleFine = new real[cellsFine];
	real* zScaleFine = new real[cellsFine];
	real* etaScaleFine = new real[cellsFine];

	real qMax = 0;
	real zMax = 0;
	real etaMax = 0;

	// finest scale coefficients and maximums
	for (i = 0; i < cellsFine; i++)
	{
		qScaleFine[i] = (qIntFine[i] + qIntFine[i + 1]) / 2;
		hScaleFine[i] = (hIntFine[i] + hIntFine[i + 1]) / 2;
		zScaleFine[i] = (zIntFine[i] + zIntFine[i + 1]) / 2;
		etaScaleFine[i] = zScaleFine[i] + hScaleFine[i];

		qMax = max(qMax, abs(qScaleFine[i]));
		zMax = max(zMax, abs(zScaleFine[i]));
		etaMax = max(etaMax, abs(etaScaleFine[i]));
	}

	int scaleCoeffsPerCell = myPowInt(2, solverParameters.L + 1) - 1;
	int totalScaleCoeffs = scaleCoeffsPerCell * simulationParameters.cells;

	real* qScaleFlattened = new real[totalScaleCoeffs];
	real* etaScaleFlattened = new real[totalScaleCoeffs];
	real* zScaleFlattened = new real[totalScaleCoeffs];
	
	real* qScaleCoarse = new real[simulationParameters.cells];
	real* etaScaleCoarse = new real[simulationParameters.cells];
	real* zScaleCoarse = new real[simulationParameters.cells];

	int detailsPerCell = myPowInt(2, solverParameters.L) - 1;
	int totalDetails = detailsPerCell * simulationParameters.cells;

	real* qDetailsFlattened = new real[totalDetails];
	real* etaDetailsFlattened = new real[totalDetails];
	real* zDetailsFlattened = new real[totalDetails];

	// START ENCODING SCALE AND DETAIL COEFFICIENTS //

	for (int c = 0; c < simulationParameters.cells; c++)
	{
		int scaleStep = c * scaleCoeffsPerCell;
		int detailStep = c * detailsPerCell;

		// loading the fine scale data
		for (int k = 0; k < myPowInt(2, solverParameters.L); k++)
		{
			qScaleFlattened[scaleStep + myPowInt(2, solverParameters.L) - 1 + k] = qScaleFine[c * myPowInt(2, solverParameters.L) + k];
			etaScaleFlattened[scaleStep + myPowInt(2, solverParameters.L) - 1 + k] = etaScaleFine[c * myPowInt(2, solverParameters.L) + k];
			zScaleFlattened[scaleStep + myPowInt(2, solverParameters.L) - 1 + k] = zScaleFine[c * myPowInt(2, solverParameters.L) + k];
		}

		for (int n = solverParameters.L - 1; n >= 0; n--)
		{
			int minLevIdx = myPowInt(2, n) - 1;
			int maxLevIdx = myPowInt(2, n + 1) - 2;
			int kHigher = maxLevIdx + 1; // index of child element

			for (int k = minLevIdx; k <= maxLevIdx; k++)
			{
				// get child elements
				real q1 = qScaleFlattened[scaleStep + kHigher];
				real eta1 = etaScaleFlattened[scaleStep + kHigher];
				real z1 = zScaleFlattened[scaleStep + kHigher];
				
				real q2 = qScaleFlattened[scaleStep + kHigher + 1];
				real eta2 = etaScaleFlattened[scaleStep + kHigher + 1];
				real z2 = zScaleFlattened[scaleStep + kHigher + 1];

				qScaleFlattened[scaleStep + k] = encodeScale(q1, q2);
				etaScaleFlattened[scaleStep + k] = encodeScale(eta1, eta2);
				zScaleFlattened[scaleStep + k] = encodeScale(z1, z2);

				qDetailsFlattened[detailStep + k] = encodeDetail(q1, q2);
				etaDetailsFlattened[detailStep + k] = encodeDetail(eta1, eta2);
				zDetailsFlattened[detailStep + k] = encodeDetail(z1, z2);

				kHigher += 2; // increment by two to step along both of the child elements
			}
		}

		qScaleCoarse[c] = qScaleFlattened[scaleStep];
		etaScaleCoarse[c] = etaScaleFlattened[scaleStep];
		zScaleCoarse[c] = zScaleFlattened[scaleStep];
	}

	// END ENCODING SCALE AND DETAIL COEFFICIENTS //

	real timeNow = 0;
	real dt = C(1e-4);

	real* normalisedDetails = new real[totalDetails];

	while (timeNow < simulationParameters.simulationTime)
	{
		timeNow += dt;

		if (timeNow - simulationParameters.simulationTime > 0)
		{
			timeNow -= dt;
			dt = simulationParameters.simulationTime - timeNow;
			timeNow += dt;
		}

		for (i = 0; i < totalDetails; i++)
		{
			qMax != 0 ? qDetailsFlattened[i] /= qMax : qDetailsFlattened[i];
			etaMax != 0 ? etaDetailsFlattened[i] /= etaMax : etaDetailsFlattened[i];
			zMax != 0 ? zDetailsFlattened[i] /= zMax : zDetailsFlattened[i];

			real a = max(qDetailsFlattened[i], etaDetailsFlattened[i]);
			normalisedDetails[i] = max(a, zDetailsFlattened[i]);
		}

		bool* significantDetails = new bool[totalDetails]();

		// START PREDICTION //

		for (int c = 0; c < simulationParameters.cells; c++)
		{
			int scaleStep = c * scaleCoeffsPerCell;
			int detailStep = c * detailsPerCell;

			for (int n = 0; n < solverParameters.L; n++)
			{
				int minLevIdx = myPowInt(2, n) - 1;
				int maxLevIdx = myPowInt(2, n + 1) - 2;

				for (int k = minLevIdx; k <= maxLevIdx; k++)
				{
					real epsilonLocal = solverParameters.epsilon * pow(C(2.0), n - solverParameters.L);

					if (normalisedDetails[detailStep + k] > epsilonLocal)
					{
						significantDetails[detailStep + k] = true;

						if (c + 1 < simulationParameters.cells && k == maxLevIdx) // if it's not the last cell and at the right-most edge
						{
							significantDetails[(c + 1) * detailsPerCell + minLevIdx] = true; // make the subelement to the right cell also significant
						}

						if (c > 0 && k == minLevIdx) // if it's not the first cell and at the left-most edge
						{
							significantDetails[(c - 1) * detailsPerCell + maxLevIdx] = true; // the subelement on the left cell is also significant
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
				int k = myPowInt(2, n - 1) - 1;
				int kLower = myPowInt(2, n - 2) - 1;
				int maxLevIdx = myPowInt(2, n) - 2;

				while (k <= maxLevIdx)
				{
					if (significantDetails[detailStep + k] || significantDetails[detailStep + k + 1])
					{
						significantDetails[detailStep + kLower] = true;
					}

					k += 2; // step along two child subelements
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
				int minLevIdx = myPowInt(2, n) - 1;
				int maxLevIdx = myPowInt(2, n + 1) - 2;
				int kHigher = maxLevIdx + 1; // index of child element

				for (int k = minLevIdx; k <= maxLevIdx; k++)
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

		delete[] qScaleFlattened;
		delete[] etaScaleFlattened;
		delete[] zScaleFlattened;

		qScaleFlattened = new real[totalScaleCoeffs];
		etaScaleFlattened = new real[totalScaleCoeffs];
		zScaleFlattened = new real[totalScaleCoeffs];

		real* dxFlattened = new real[totalScaleCoeffs];
		real* xFlattened = new real[totalScaleCoeffs];
		int* levelIndicesFlattened = new int[totalScaleCoeffs];

		vector<real> qAdaptive;
		vector<real> etaAdaptive;
		vector<real> zAdaptive;
		vector<real> dxAdaptive;
		vector<real> xAdaptive;
		vector<int> levelIndicesAdaptive;
		vector<int> activeIndices;

		// initialising with NAN to later check for undecoded scale coeff and avoid assembly
		for (int i = 0; i < totalScaleCoeffs; i++)
		{
			qScaleFlattened[i] = NAN;
		}

		real* xIntCoarse = new real[simulationParameters.cells + 1];

		// coarse mesh
		for (i = 0; i < simulationParameters.cells + 1; i++)
		{
			xIntCoarse[i] = simulationParameters.xmin + dxCoarse * i;
		}

		for (int c = 0; c < simulationParameters.cells; c++)
		{
			int scaleStep = c * scaleCoeffsPerCell;

			qScaleFlattened[scaleStep] = qScaleCoarse[c];
			etaScaleFlattened[scaleStep] = etaScaleCoarse[c];
			zScaleFlattened[scaleStep] = zScaleCoarse[c];
			xFlattened[scaleStep] = (xIntCoarse[c] + xIntCoarse[c + 1]) / 2;
			dxFlattened[scaleStep] = dxCoarse;
			levelIndicesFlattened[scaleStep] = 0;
		}

		// START DECODING //

		for (int c = 0; c < simulationParameters.cells; c++)
		{
			int scaleStep = c * scaleCoeffsPerCell;
			int detailStep = c * detailsPerCell;

			for (int n = 0; n < solverParameters.L; n++)
			{
				int minLevIdx = myPowInt(2, n) - 1;
				int maxLevIdx = myPowInt(2, n + 1) - 2;
				int kHigher = maxLevIdx + 1; // index of child element

				for (int k = minLevIdx; k <= maxLevIdx; k++)
				{
					if (!significantDetails[detailStep + k] && !isnan(qScaleFlattened[scaleStep + k]))
					{
						real q = qScaleFlattened[scaleStep + k];
						real eta = etaScaleFlattened[scaleStep + k];
						real z = zScaleFlattened[scaleStep + k];
						real dxLocal = dxFlattened[scaleStep + k];
						real x = xFlattened[scaleStep + k];
						int levelIndex = levelIndicesFlattened[scaleStep + k];

						qAdaptive.push_back(q);
						etaAdaptive.push_back(eta);
						zAdaptive.push_back(z);
						dxAdaptive.push_back(dxLocal);
						xAdaptive.push_back(x);
						levelIndicesAdaptive.push_back(levelIndex);
						activeIndices.push_back(scaleStep + k);
					}
					else if (significantDetails[detailStep + k])
					{
						real q1 = decode1(qScaleFlattened[scaleStep + k], qDetailsFlattened[detailStep + k]);
						real eta1 = decode1(etaScaleFlattened[scaleStep + k], etaDetailsFlattened[detailStep + k]);
						real z1 = decode1(zScaleFlattened[scaleStep + k], zDetailsFlattened[detailStep + k]);

						real q2 = decode2(qScaleFlattened[scaleStep + k], qDetailsFlattened[detailStep + k]);
						real eta2 = decode2(etaScaleFlattened[scaleStep + k], etaDetailsFlattened[detailStep + k]);
						real z2 = decode2(zScaleFlattened[scaleStep + k], zDetailsFlattened[detailStep + k]);

						qScaleFlattened[scaleStep + kHigher] = q1;
						etaScaleFlattened[scaleStep + kHigher] = eta1;
						zScaleFlattened[scaleStep + kHigher] = z1;

						qScaleFlattened[scaleStep + kHigher + 1] = q2;
						etaScaleFlattened[scaleStep + kHigher + 1] = eta2;
						zScaleFlattened[scaleStep + kHigher + 1] = z2;

						real dxHigher = dxFlattened[scaleStep + k] / 2;
						dxFlattened[scaleStep + kHigher] = dxHigher;
						dxFlattened[scaleStep + kHigher + 1] = dxHigher;

						real x1 = xFlattened[scaleStep + k] - dxFlattened[scaleStep + kHigher] / 2;
						real x2 = xFlattened[scaleStep + k] + dxFlattened[scaleStep + kHigher] / 2;
						xFlattened[scaleStep + kHigher] = x1;
						xFlattened[scaleStep + kHigher + 1] = x2;

						int nHigher = levelIndicesFlattened[scaleStep + k] + 1;
						levelIndicesFlattened[scaleStep + kHigher] = nHigher;
						levelIndicesFlattened[scaleStep + kHigher + 1] = nHigher;

						if (levelIndicesFlattened[scaleStep + k] + 1 == solverParameters.L)
						{
							qAdaptive.push_back(q1);
							qAdaptive.push_back(q2);

							etaAdaptive.push_back(eta1);
							etaAdaptive.push_back(eta2);

							zAdaptive.push_back(z1);
							zAdaptive.push_back(z2);

							dxAdaptive.push_back(dxHigher);
							dxAdaptive.push_back(dxHigher);

							xAdaptive.push_back(x1);
							xAdaptive.push_back(x2);

							levelIndicesAdaptive.push_back(nHigher);
							levelIndicesAdaptive.push_back(nHigher);

							activeIndices.push_back(scaleStep + kHigher);
							activeIndices.push_back(scaleStep + kHigher + 1);
						}
					}

					kHigher += 2;
				}
			}
		}

		// END DECODING //


		// START SORTING //

		vector<real> xSorted = xAdaptive;
		vector<pair<real, real>> qZipped;
		zip(qAdaptive, xSorted, qZipped);

		// sort the vector of pairs
		sort(begin(qZipped), end(qZipped),
			[&](const auto& a, const auto& b)
			{
				return a.second < b.second;
			});

		// write the sorted pairs back to the original vectors
		unzip(qZipped, qAdaptive, xSorted);

		xSorted = xAdaptive;
		vector<pair<real, real>> etaZipped;
		zip(etaAdaptive, xSorted, etaZipped);

		sort(begin(etaZipped), end(etaZipped),
			[&](const auto& a, const auto& b)
			{
				return a.second < b.second;
			});

		unzip(etaZipped, etaAdaptive, xSorted);

		xSorted = xAdaptive;
		vector<pair<real, real>> zZipped;
		zip(zAdaptive, xSorted, zZipped);

		sort(begin(zZipped), end(zZipped),
			[&](const auto& a, const auto& b)
			{
				return a.second < b.second;
			});

		unzip(zZipped, zAdaptive, xSorted);

		xSorted = xAdaptive;
		vector<pair<real, real>> dxZipped;
		zip(dxAdaptive, xSorted, dxZipped);

		sort(begin(dxZipped), end(dxZipped),
			[&](const auto& a, const auto& b)
			{
				return a.second < b.second;
			});

		unzip(dxZipped, dxAdaptive, xSorted);

		xSorted = xAdaptive;
		vector<pair<int, real>> levelIndicesZipped;
		zip(levelIndicesAdaptive, xSorted, levelIndicesZipped);

		sort(begin(levelIndicesZipped), end(levelIndicesZipped),
			[&](const auto& a, const auto& b)
			{
				return a.second < b.second;
			});

		unzip(levelIndicesZipped, levelIndicesAdaptive, xSorted);

		xSorted = xAdaptive;
		vector<pair<int, real>> activeIndicesZipped;
		zip(activeIndices, xSorted, activeIndicesZipped);

		sort(begin(activeIndicesZipped), end(activeIndicesZipped),
			[&](const auto& a, const auto& b)
			{
				return a.second < b.second;
			});

		unzip(activeIndicesZipped, activeIndices, xSorted);

		// END SORTING //


		int assembledSolutionLength = activeIndices.size();

		// allocate buffers for flow modes with ghost BCs
		real* qWithBC = new real[assembledSolutionLength + 2];
		real* hWithBC = new real[assembledSolutionLength + 2];
		real* zWithBC = new real[assembledSolutionLength + 2];
		real* dxLocalWithBC = new real[assembledSolutionLength + 2];

		// load from vector 
		for (i = 0; i < assembledSolutionLength; i++)
		{
			qWithBC[i + 1] = qAdaptive[i];
			hWithBC[i + 1] = etaAdaptive[i] - zAdaptive[i];
			zWithBC[i + 1] = zAdaptive[i];
			dxLocalWithBC[i + 1] = xSorted[i];
		}

		// allocate true/false buffer for dry cells
		bool* dry = new bool[assembledSolutionLength + 2]();

		real* etaTemp = new real[assembledSolutionLength + 2];

		// allocating buffers for eastern and western interface values
		real* qEast = new real[assembledSolutionLength + 1];
		real* hEast = new real[assembledSolutionLength + 1];
		real* etaEast = new real[assembledSolutionLength + 1];

		real* qWest = new real[assembledSolutionLength + 1];
		real* hWest = new real[assembledSolutionLength + 1];
		real* etaWest = new real[assembledSolutionLength + 1];

		// allocating buffers for positivity preserving nodes
		real* qEastStar = new real[assembledSolutionLength + 1];
		real* hEastStar = new real[assembledSolutionLength + 1];

		real* qWestStar = new real[assembledSolutionLength + 1];
		real* hWestStar = new real[assembledSolutionLength + 1];

		real* zStarIntermediate = new real[assembledSolutionLength + 1];
		real* zStar = new real[assembledSolutionLength + 1];

		real* uWest = new real[assembledSolutionLength + 1];
		real* uEast = new real[assembledSolutionLength + 1];

		real* delta = new real[assembledSolutionLength + 1];

		// allocating buffers for numerical fluxes from HLL solver
		real* massFlux = new real[assembledSolutionLength + 1];
		real* momentumFlux = new real[assembledSolutionLength + 1];

		// allocating buffers for positivity preserving MODES
		real* hBar = new real[assembledSolutionLength];
		real* zBar = new real[assembledSolutionLength];

		// adding ghost boundary conditions
		real qUp = bcs.reflectUp * qWithBC[1];
		if (bcs.qxImposedUp > 0)
		{
			qUp = bcs.qxImposedUp;
		}

		real qDown = bcs.reflectDown * qWithBC[assembledSolutionLength];
		if (bcs.qxImposedDown > 0)
		{
			qDown = bcs.qxImposedDown;
		}

		qWithBC[0] = qUp;
		qWithBC[assembledSolutionLength + 1] = qDown; // recall there are cells + 2 elements inc BCs

		real hUp = hWithBC[1];
		if (bcs.hImposedUp > 0)
		{
			hUp = bcs.hImposedUp;
		}

		real hDown = hWithBC[assembledSolutionLength];
		if (bcs.hImposedDown > 0)
		{
			hDown = bcs.hImposedDown;
		}

		hWithBC[0] = hUp;
		hWithBC[assembledSolutionLength + 1] = hDown;

		real zUp = zWithBC[1];
		real zDown = zWithBC[assembledSolutionLength];

		zWithBC[0] = zUp;
		zWithBC[assembledSolutionLength + 1] = zDown;

		dxLocalWithBC[0] = dxLocalWithBC[1];
		dxLocalWithBC[assembledSolutionLength + 1] = dxLocalWithBC[assembledSolutionLength];

		// extract upwind and downwind modes
		real hWestUpwind = hWithBC[0];
		real hEastDownwind = hWithBC[assembledSolutionLength + 1];

		if (simulationParameters.manning > 0)
		{
			for (i = 1; i < assembledSolutionLength + 1; i++)
			{
				qWithBC[i] += frictionImplicit(simulationParameters, solverParameters, dt, hWithBC[i], qWithBC[i]);
			}
		}

		// initialising dry vs wet cells, ignore ghost cells
		for (i = 1; i < assembledSolutionLength + 1; i++)
		{
			real hLocal = hWithBC[i];
			real hBackward = hWithBC[i - 1];
			real hForward = hWithBC[i + 1];

			real hMax = max(hLocal, hBackward);
			hMax = max(hForward, hMax);

			if (hMax <= solverParameters.tolDry)
			{
				dry[i] = true;
			}
		}

		for (i = 0; i < assembledSolutionLength + 2; i++)
		{
			etaTemp[i] = hWithBC[i] + zWithBC[i];
		}

		// initialising interface values
		for (i = 0; i < assembledSolutionLength + 1; i++)
		{
			qEast[i] = qWithBC[i + 1];
			hEast[i] = hWithBC[i + 1];
			etaEast[i] = etaTemp[i + 1];

			qWest[i] = qWithBC[i];
			hWest[i] = hWithBC[i];
			etaWest[i] = etaTemp[i];
		}

		// correcting downwind and upwind eta values
		etaEast[assembledSolutionLength] = etaTemp[assembledSolutionLength] - hWithBC[assembledSolutionLength] + hEastDownwind;
		etaWest[0] = etaTemp[1] - hWithBC[1] + hWestUpwind;

		for (int i = 0; i < assembledSolutionLength + 1; i++)
		{
			// initialising velocity interface values
			uWest[i] = uFaceValues(solverParameters, qWest[i], hWest[i]);
			uEast[i] = uFaceValues(solverParameters, qEast[i], hEast[i]);

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
		fluxHLL(assembledSolutionLength, solverParameters, hWestStar, hEastStar, qWestStar, qEastStar, uWest, uEast, massFlux, momentumFlux);

		for (int i = 0; i < assembledSolutionLength; i++)
		{
			// essentially 0th order projection but taking into account east/west locality
			hBar[i] = (hEastStar[i] + hWestStar[i + 1]) / 2;

			// 1st order projection
			zBar[i] = (zStar[i + 1] - zStar[i]) / (2 * sqrt(C(3.0)));
		}

		// FV1 operator increment, skip ghosts cells
		for (i = 1; i < assembledSolutionLength + 1; i++)
		{
			// skip increment in dry cells
			if (dry[i])
			{
				continue;
			}
			else
			{
				real massIncrement = -(1 / dxLocalWithBC[i]) * (massFlux[i] - massFlux[i - 1]);
				real momentumIncrement = -(1 / dxLocalWithBC[i]) * (momentumFlux[i] - momentumFlux[i - 1] + 2 * sqrt(C(3.0)) * solverParameters.g * hBar[i - 1] * zBar[i - 1]);

				hWithBC[i] += dt * massIncrement;
				qWithBC[i] += dt * momentumIncrement;
			}

			if (hWithBC[i] <= solverParameters.tolDry)
			{
				qWithBC[i] = 0;
			}
		}

		// CFL time step adjustment
		dt = 1e9;

		for (i = 1; i < assembledSolutionLength + 1; i++)
		{
			if (hWithBC[i] <= solverParameters.tolDry)
			{
				continue;
			}
			else
			{
				real u = qWithBC[i] / hWithBC[i];
				real dtCFL = solverParameters.CFL * dxLocalWithBC[i] / (abs(u) + sqrt(solverParameters.g * hWithBC[i]));
				dt = min(dt, dtCFL);
			}
		}

		for (i = 1; i < assembledSolutionLength + 1; i++)
		{
			printf("%.2f, ", hWithBC[i] + zWithBC[i]);
		}
		printf("\n");

		delete[] qScaleFlattened;
		delete[] etaScaleFlattened;
		delete[] zScaleFlattened;

		qScaleFlattened = new real[totalScaleCoeffs];
		etaScaleFlattened = new real[totalScaleCoeffs];
		zScaleFlattened = new real[totalScaleCoeffs];

		delete[] qDetailsFlattened;
		delete[] etaDetailsFlattened;
		delete[] zDetailsFlattened;

		// initialise zeroed details for thresholding
		qDetailsFlattened = new real[totalDetails]();
		etaDetailsFlattened = new real[totalDetails]();
		zDetailsFlattened = new real[totalDetails]();

		qMax = 0;
		etaMax = 0;
		zMax = 0;

		for (i = 0; i < assembledSolutionLength; i++)
		{
			real q = qWithBC[i + 1];
			real eta = hWithBC[i + 1] + zWithBC[i + 1];
			real z = zWithBC[i + 1];

			qScaleFlattened[activeIndices[i]] = q;
			etaScaleFlattened[activeIndices[i]] = eta;
			zScaleFlattened[activeIndices[i]] = z;

			qMax = max(q, qMax);
			etaMax = max(eta, etaMax);
			zMax = max(z, zMax);
		}

		
		// BEGIN ENCODING //

		for (int c = 0; c < simulationParameters.cells; c++)
		{
			int scaleStep = c * scaleCoeffsPerCell;
			int detailStep = c * detailsPerCell;

			for (int n = solverParameters.L - 1; n >= 0; n--)
			{
				int minLevIdx = myPowInt(2, n) - 1;
				int maxLevIdx = myPowInt(2, n + 1) - 2;
				int kHigher = maxLevIdx + 1;

				for (int k = minLevIdx; k <= maxLevIdx; k++)
				{
					if (significantDetails[detailStep + k])
					{
						real q1 = qScaleFlattened[scaleStep + kHigher];
						real eta1 = etaScaleFlattened[scaleStep + kHigher];
						real z1 = zScaleFlattened[scaleStep + kHigher];

						real q2 = qScaleFlattened[scaleStep + kHigher + 1];
						real eta2 = etaScaleFlattened[scaleStep + kHigher + 1];
						real z2 = zScaleFlattened[scaleStep + kHigher + 1];

						qScaleFlattened[scaleStep + k] = encodeScale(q1, q2);
						etaScaleFlattened[scaleStep + k] = encodeScale(eta1, eta2);
						zScaleFlattened[scaleStep + k] = encodeScale(z1, z2);

						qDetailsFlattened[detailStep + k] = encodeDetail(q1, q2);
						etaDetailsFlattened[detailStep + k] = encodeDetail(eta1, eta2);
						zDetailsFlattened[detailStep + k] = encodeDetail(z1, z2);
					}

					kHigher += 2;
				}
			}

			qScaleCoarse[c] = qScaleFlattened[scaleStep];
			etaScaleCoarse[c] = etaScaleFlattened[scaleStep];
			zScaleCoarse[c] = zScaleFlattened[scaleStep];
		}

		// END ENCODING //

		//for (int c = 0; c < simulationParameters.cells; c++)
		//{
		//	int detailStep = c * detailsPerCell;

		//	for (int n = solverParameters.L - 1; n >= 0; n--)
		//	{
		//		int minLevIdx = myPowInt(2, n) - 1;
		//		int maxLevIdx = myPowInt(2, n + 1) - 2;
		//		int kHigher = maxLevIdx + 1; // index of child element

		//		for (int k = minLevIdx; k <= maxLevIdx; k++)
		//		{
		//			printf(" %d |", significantDetails[detailStep + k] ? 1 : 0);
		//		}

		//		printf("\n");
		//	}

		//	printf("\n");
		//}

		printf("%f s\n", timeNow);

		delete[] significantDetails;

		delete[] dxFlattened;
		delete[] xFlattened;
		delete[] levelIndicesFlattened;

		delete[] xIntCoarse;

		delete[] qWithBC;
		delete[] hWithBC;
		delete[] zWithBC;
		delete[] dxLocalWithBC;

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

	delete[] qScaleFlattened;
	delete[] etaScaleFlattened;
	delete[] zScaleFlattened;

	delete[] qScaleCoarse;
	delete[] etaScaleCoarse;
	delete[] zScaleCoarse;

	delete[] qDetailsFlattened;
	delete[] etaDetailsFlattened;
	delete[] zDetailsFlattened;

	delete[] normalisedDetails;

	return 0;
}

template <typename A, typename B> // typename refers to any generic type
void zip(const vector<A>& a, const vector<B>& b, vector<pair<A, B>>& zipped)
{
	for (size_t i = 0; i < a.size(); ++i)
	{
		zipped.push_back(make_pair(a[i], b[i]));
	}
}

template <typename A, typename B>
void unzip(const vector<pair<A, B>>& zipped, vector<A>& a, vector<B>& b)
{
	for (size_t i = 0; i < a.size(); i++)
	{
		a[i] = zipped[i].first;
		b[i] = zipped[i].second;
	}
}