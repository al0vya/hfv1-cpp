#include "writeSolutionToFile.h"

void writeSolutionToFile(AssembledSolution assembledSolution)
{
	std::ofstream test;

	test.open("debug.csv");

	test.precision(17);

	test << "x,q,z,eta" << std::endl;

	for (int i = 0; i < assembledSolution.length; i++)
	{
		test << assembledSolution.x[i + 1] << "," << assembledSolution.qWithBC[i + 1] << "," << assembledSolution.zWithBC[i + 1] << "," << 
			std::max(assembledSolution.hWithBC[i + 1] + assembledSolution.zWithBC[i + 1], assembledSolution.hWithBC[i + 1]) << std::endl;
	}

	test.close();
}