#pragma once
#include <cmath>
#include <algorithm>

using namespace std;

#include "real.h"
#include "SolverParameters.h"
#include "AssembledSolution.h"

void fluxHLL(AssembledSolution assembledSolution, SolverParameters solverParameters, real* hWestStar, real* hEastStar, real* qWestStar, real* qEastStar, real* massFlux, real* momentumFlux);