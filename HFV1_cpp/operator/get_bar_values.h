#pragma once

#include <algorithm>

#include "../classes/StarValues.h"
#include "../classes/AssembledSolution.h"
#include "../classes/BarValues.h"

void get_bar_values
(
	StarValues&        star_vals, 
	AssembledSolution& assem_sol, 
	BarValues&         bar_vals
);