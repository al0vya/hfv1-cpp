#include "malloc_scale_coefficients.h"

void malloc_scale_coefficients
(
	FlattenedScaleCoeffs& scale_coeffs,
	int&                  num_scale_coeffs
)
{
	scale_coeffs.q   = new real[num_scale_coeffs]();
	scale_coeffs.eta = new real[num_scale_coeffs]();
	scale_coeffs.z   = new real[num_scale_coeffs]();
}