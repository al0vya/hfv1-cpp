#pragma once

#include "real.h"

typedef struct FlattenedScaleCoeffs
{
	real* q;
	real* eta;
	real* z;

	bool is_copy = false;

	FlattenedScaleCoeffs(const int& num_scale_coeffs)
	{
		q   = new real[num_scale_coeffs]();
	    eta = new real[num_scale_coeffs]();
	    z   = new real[num_scale_coeffs]();
	}

	FlattenedScaleCoeffs(const FlattenedScaleCoeffs& original) { *this = original; is_copy = true; }

	~FlattenedScaleCoeffs()
	{
		if (!is_copy)
		{
			delete[] q;
			delete[] eta;
			delete[] z;
		}
	}

} FlattenedScaleCoeffs;