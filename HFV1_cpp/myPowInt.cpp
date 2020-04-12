#include "myPowInt.h"

int myPowInt(int base, int exponent)
{
	int a = base;

	if (exponent == 0)
	{
		return 1;
	}
	else
	{
		for (int i = 1; i < exponent; i++)
		{
			a *= base;
		}
	}

	return a;
}