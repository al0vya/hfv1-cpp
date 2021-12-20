#pragma once

#include "AssembledSolution.h"
#include "SimulationParameters.h"
#include "FaceValues.h"

void get_face_values
(
	AssembledSolution& assem_sol, 
	FaceValues&        face_vals, 
	real*&             eta_temp
);