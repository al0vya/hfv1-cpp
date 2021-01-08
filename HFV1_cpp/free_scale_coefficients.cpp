#include "free_scale_coefficients.h"

void free_scale_coefficients(FlattenedScaleCoeffs& scale_coeffs)
{
	delete[] scale_coeffs.q;
	delete[] scale_coeffs.eta;
	delete[] scale_coeffs.z;
}