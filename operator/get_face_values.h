#pragma once

#include "../classes/AssembledSolution.h"
#include "../classes/SimulationParameters.h"
#include "../classes/FaceValues.h"

void get_face_values
(
	AssembledSolution& assem_sol, 
	FaceValues&        face_vals, 
	real*&             eta_temp
);