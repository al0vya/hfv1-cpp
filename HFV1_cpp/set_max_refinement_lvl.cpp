#include "set_max_refinement_lvl.h"

int set_max_refinement_lvl()
{
	std::cout << "Please select a maximum refinement level.\n";

	int user_input_max_refinement_level;

	std::cin >> user_input_max_refinement_level;

	if (!std::cin || user_input_max_refinement_level < 1)
	{
		std::cout << "Error: please rerun and enter a integer value. Exiting program.\n";

		exit(0);
	}

	return user_input_max_refinement_level;
}