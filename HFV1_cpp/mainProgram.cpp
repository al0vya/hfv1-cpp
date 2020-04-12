#include <stdlib.h>
#include <stdio.h>
#include <algorithm>
#include <cmath>

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

using namespace std;

int main()
{
	// quintessential for-loop index
	int i;

	int step = 0;

	SimulationParameters simulationParameters;
	simulationParameters.cells = 2;
	simulationParameters.xmin = 0;
	simulationParameters.xmax = 50;
	simulationParameters.simulationTime = C(20.);
	simulationParameters.manning = C(0.02);

	SolverParameters solverParameters;
	solverParameters.CFL = C(0.33);
	solverParameters.tolDry = C(1e-3);
	solverParameters.g = C(9.80665);
	solverParameters.L = 7;

	BoundaryConditions bcs;
	bcs.hl = C(6.0);
	bcs.hr = C(1.0);

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
	real* qScaleFlattened = new real[simulationParameters.cells * scaleCoeffsPerCell];
	real* etaScaleFlattened = new real[simulationParameters.cells * scaleCoeffsPerCell];
	real* zScaleFlattened = new real[simulationParameters.cells * scaleCoeffsPerCell];
	
	real* qScaleCoarse = new real[simulationParameters.cells];
	real* etaScaleCoarse = new real[simulationParameters.cells];
	real* zScaleCoarse = new real[simulationParameters.cells];

	int detailsPerCell = myPowInt(2, solverParameters.L) - 1;
	real* qDetailsFlattened = new real[simulationParameters.cells * detailsPerCell];
	real* etaDetailsFlattened = new real[simulationParameters.cells * detailsPerCell];
	real* zDetailsFlattened = new real[simulationParameters.cells * detailsPerCell];

	// START ENCODING SCALE AND DETAIL COEFFICIENTS //

	for (int c = 0; c < simulationParameters.cells; c++)
	{
		for (int k = 0; k < myPowInt(2, solverParameters.L); k++)
		{
			qScaleFlattened[c * scaleCoeffsPerCell + myPowInt(2, solverParameters.L) - 1 + k] = qScaleFine[c * myPowInt(2, solverParameters.L) + k];
			etaScaleFlattened[c * scaleCoeffsPerCell + myPowInt(2, solverParameters.L) - 1 + k] = etaScaleFine[c * myPowInt(2, solverParameters.L) + k];
			zScaleFlattened[c * scaleCoeffsPerCell + myPowInt(2, solverParameters.L) - 1 + k] = zScaleFine[c * myPowInt(2, solverParameters.L) + k];
		}

		for (int n = solverParameters.L - 1; n >= 0; n--)
		{
			int minLevIdx = myPowInt(2, n) - 1;
			int maxLevIdx = myPowInt(2, n + 1) - 2;
			int kHigher = myPowInt(2, n + 1) - 1; // index of child element

			for (int k = minLevIdx; k <= maxLevIdx; k++)
			{
				// get child elements
				real q1 = qScaleFlattened[c * scaleCoeffsPerCell + kHigher];
				real eta1 = etaScaleFlattened[c * scaleCoeffsPerCell + kHigher];
				real z1 = zScaleFlattened[c * scaleCoeffsPerCell + kHigher];
				
				real q2 = qScaleFlattened[c * scaleCoeffsPerCell + kHigher + 1];
				real eta2 = etaScaleFlattened[c * scaleCoeffsPerCell + kHigher + 1];
				real z2 = zScaleFlattened[c * scaleCoeffsPerCell + kHigher + 1];

				qScaleFlattened[c * scaleCoeffsPerCell + k] = encodeScale(q1, q2);
				etaScaleFlattened[c * scaleCoeffsPerCell + k] = encodeScale(eta1, eta2);
				zScaleFlattened[c * scaleCoeffsPerCell + k] = encodeScale(z1, z2);

				qDetailsFlattened[c * detailsPerCell + k] = encodeDetail(q1, q2);
				etaDetailsFlattened[c * detailsPerCell + k] = encodeDetail(eta1, eta2);
				zDetailsFlattened[c * detailsPerCell + k] = encodeDetail(z1, z2);

				kHigher += 2; // increment by two to step along both of the child elements
			}
		}

		qScaleCoarse[c] = qScaleFlattened[c * scaleCoeffsPerCell];
		etaScaleCoarse[c] = etaScaleFlattened[c * scaleCoeffsPerCell];
		zScaleCoarse[c] = zScaleFlattened[c * scaleCoeffsPerCell];
	}

	// END ENCODING SCALE AND DETAIL COEFFICIENTS //

	real* normalisedDetails = new real[simulationParameters.cells * detailsPerCell];

	for (i = 0; i < simulationParameters.cells * detailsPerCell; i++)
	{
		qMax != 0 ? qDetailsFlattened[i] /= qMax : qDetailsFlattened[i];
		etaMax != 0 ? etaDetailsFlattened[i] /= etaMax : etaDetailsFlattened[i];
		zMax != 0 ? zDetailsFlattened[i] /= zMax : zDetailsFlattened[i];

		real a = max(qDetailsFlattened[i], etaDetailsFlattened[i]);
		normalisedDetails[i] = max(a, zDetailsFlattened[i]);
	}

	bool* significantDetails = new bool[simulationParameters.cells * detailsPerCell]();

	// START PREDICTION //
	
	for (int c = 0; c < simulationParameters.cells; c++)
	{
		for (int n = 0; n < solverParameters.L; n++)
		{
			int minLevIdx = myPowInt(2, n) - 1;
			int maxLevIdx = myPowInt(2, n + 1) - 2;
			int kHigher = myPowInt(2, n + 1) - 1; // index of child element

			for (int k = minLevIdx; k <= maxLevIdx; k++)
			{
				if (normalisedDetails[c * detailsPerCell + k] > solverParameters.epsilon * pow(C(2.0), n - solverParameters.L))
				{
					significantDetails[c * detailsPerCell + k] = true;

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
		for (int n = solverParameters.L; n > 1; n--)
		{
			int k = myPowInt(2, n - 1) - 1;
			int kLower = myPowInt(2, n - 2) - 1;

			while (k < myPowInt(2, n) - 2)
			{
				if (significantDetails[c * detailsPerCell + k] || significantDetails[c * detailsPerCell + k + 1])
				{
					significantDetails[c * detailsPerCell + kLower] = true;
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
		for (int n = 0; n < solverParameters.L; n++)
		{
			int minLevIdx = myPowInt(2, n) - 1;
			int maxLevIdx = myPowInt(2, n + 1) - 2;
			int kHigher = myPowInt(2, n + 1) - 1; // index of child element

			for (int k = minLevIdx; k <= maxLevIdx; k++)
			{
				if (significantDetails[c * detailsPerCell + k])
				{
					real mBar = 1.5;
					real epsilonLocal = solverParameters.epsilon * pow(C(2.0), n - solverParameters.L);

					if (normalisedDetails[c * detailsPerCell + k] >= epsilonLocal * pow(2, mBar + 1) && n + 1 != solverParameters.L)
					{
						// if extra signficant child elements marked as active
						significantDetails[c * detailsPerCell + kHigher] = true;
						significantDetails[c * detailsPerCell + kHigher + 1] = true;
					}
				}

				kHigher += 2;
			}
		}
	}

	// END EXTRA SIGNIFICANCE //

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

	delete[] significantDetails;

	return 0;
}