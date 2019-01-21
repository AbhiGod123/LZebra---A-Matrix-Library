#include "Standalone.h"

float gaussianRandom() {
	float v1, v2, s;
	do {
		v1 = 2.0f * (float)(rand() % 10000) / (10000) - 1.0f;
		v2 = 2.0f * (float)(rand() % 10000) / (10000) - 1.0f;
		s = v1 * v1 + v2 * v2;
	} while (s >= 1.0f || s == 0.0f);

	s = (float)pow((-2.0f * log(s)) / s, 0.5f);

	return v1 * s;
}

float uniformFloatRandom()
{
	return (float)(rand() % 10000) / (10000);
}

int uniformIntRandom()
{
	return rand() % INT_MAX;
}