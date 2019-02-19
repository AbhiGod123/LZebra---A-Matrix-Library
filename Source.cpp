#include "LZebra.h"

//An example program that applies uniform float randoms to each element in the matrix
int main() {
	LZebra::Matrix<float> d(3, 3);

	d.randu();

	std::cout << d << '\n';

	{
		system("PAUSE");
		return 0;
	}
}