#include "free_details.h"

void free_details(FlattenedDetails& details)
{
	delete[] details.q;
	delete[] details.eta;
	delete[] details.z;
}