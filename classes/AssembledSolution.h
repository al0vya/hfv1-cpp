#pragma once

#include "real.h"

#include <cstdio>

typedef struct AssembledSolution
{
	real* q_BC;
	real* h_BC;
	real* z_BC;
	real* dx_BC;
	real* x;
	int*  activeIndices;
	int   length;
	bool  is_copy = false;

	AssembledSolution(const int& num_cells)
	:
		length(num_cells)
	{
		q_BC  = new real[num_cells + 2];
		h_BC  = new real[num_cells + 2];
		z_BC  = new real[num_cells + 2];
		dx_BC = new real[num_cells + 2];
		x     = new real[num_cells + 2];
		
		activeIndices = new int[num_cells + 2];
	}

	AssembledSolution(const AssembledSolution& original) { *this = original; is_copy = true; }

	~AssembledSolution()
	{
		if (!is_copy)
		{
			delete[] q_BC;
			delete[] h_BC;
			delete[] z_BC;
			delete[] dx_BC;
			delete[] x;
			delete[] activeIndices;
		}
	}

} AssembledSolution;