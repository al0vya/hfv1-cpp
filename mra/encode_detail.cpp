#include "encode_detail.h"

real encode_detail(real u1, real u2)
{
	//return C(0.5) * u1 - C(0.5) * u2;
	return sqrt(1 / C(2.0)) * (G0 * u1 + G1 * u2);
}