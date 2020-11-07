#include "decode1.h"

real decode1(real u1, real d1)
{
	//return u1 + d1;
	return sqrt(C(2.0)) * (H0 * u1 + G0 * d1);
}