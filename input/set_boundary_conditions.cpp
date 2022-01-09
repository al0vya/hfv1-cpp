#include "set_boundary_conditions.h"

BoundaryConditions set_boundary_conditions(int test_case)
{
	BoundaryConditions bcs;

	switch (test_case)
	{
	case 1:
		bcs.hl = C(6.0);
		bcs.hr = C(2.0);
		break;
	case 2:
		bcs.hl = C(6.0);
		bcs.hr = C(0.0);
		break;
	case 3:
		bcs.hl = C(6.0);
		bcs.hr = C(0.0);
		break;
	case 4:
		bcs.hl = C(6.0);
		bcs.hr = C(6.0);
		break;
	case 5:
		bcs.hl = C(2.0);
		bcs.hr = C(2.0);
		break;
	case 6:
		bcs.hl = C(6.0);
		bcs.hr = C(0.0);
		break;
	case 7:
		bcs.hl = C(0.75);
		bcs.hr = C(0.0);
		break;
	default:
		break;
	}

	return bcs;
}