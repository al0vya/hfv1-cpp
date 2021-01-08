#include "malloc_details.h"

void malloc_details
(
	FlattenedDetails& details, 
	int&              num_details
)
{
	details.q   = new real[num_details]();
	details.eta = new real[num_details]();
	details.z   = new real[num_details]();
}