#pragma once

#include "../classes/real.h"

typedef struct NodalValues
{
	real* q;
	real* h;
	real* z;
	real* x;
	bool is_copy = false;

	NodalValues(const int& num_interfaces)
	{
		q = new real[num_interfaces];
		h = new real[num_interfaces];
		z = new real[num_interfaces];
		x = new real[num_interfaces];
	}

	NodalValues(const NodalValues& original) { *this = original; is_copy = true; }

	~NodalValues()
	{
		if (!is_copy)
		{
			delete[] q;
			delete[] h;
			delete[] z;
			delete[] x;
		}
	}

} NodalValues;