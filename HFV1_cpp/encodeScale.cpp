#include "encodeScale.h"

real encodeScale(real u1, real u2)
{
	return C(0.5) * u1 + C(0.5) * u2;
	//return sqrt(1 / C(2.0)) * (H0 * u1 + H1 * u2);
}