#pragma once

#include "real.h"

typedef struct FlattenedDetails
{
	real* q;
	real* eta;
	real* z;

	bool is_copy = false;

	FlattenedDetails(const int& num_details)
	{
		q   = new real[num_details];
		eta = new real[num_details];
		z   = new real[num_details];
	}

	FlattenedDetails(const FlattenedDetails& original) { *this = original; is_copy = true; }

	~FlattenedDetails()
	{
		if (!is_copy)
		{
			delete[] q;
			delete[] eta;
			delete[] z;
		}
	}

} FlattenedDetails;