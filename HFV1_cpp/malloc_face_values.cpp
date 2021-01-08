#include "malloc_face_values.h"

void malloc_face_values
(
	FaceValues& face_values,
	int&        num_fine_cells
)
{
	face_values.q_east   = new real[num_fine_cells + 1];
	face_values.h_east   = new real[num_fine_cells + 1];
	face_values.eta_east = new real[num_fine_cells + 1];

	face_values.q_west   = new real[num_fine_cells + 1];
	face_values.h_west   = new real[num_fine_cells + 1];
	face_values.eta_west = new real[num_fine_cells + 1];
}