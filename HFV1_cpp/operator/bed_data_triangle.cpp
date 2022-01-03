#include "bed_data_triangle.h"

// from https://www.sciencedirect.com/science/article/pii/S0309170815002237
real bed_data_triangle
(
	const real& x
)
{
	real slope = C(0.4) / C(3.0);

	real z_int = C(0.0);

	if ( x > C(25.5) && x <= C(28.5) )
	{
		z_int = ( x - C(25.5) ) * slope;
	}
	else if ( x > C(28.5) && x <= C(31.5) )
	{
		z_int = C(0.4) - ( x - C(28.5) ) * slope;
	}

	return z_int;
}