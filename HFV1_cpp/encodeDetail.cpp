#include "encodeDetail.h"

real encodeDetail(real u1, real u2)
{
	return sqrt(1 / C(2.0)) * (G0 * u1 + G1 * u2);
}