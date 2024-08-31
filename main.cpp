#include <iostream>
#include <vector>
#include "Header files/matrix.h"




int main() {

	std::vector<Vector> matrix = {
			{
				Vector({4, -3 , 0}),
				Vector({1, 2 , 0}),
				Vector({0, 0 , 4})

			}
	};
	Matrix m(matrix);
	std::cout<<m.rref().toString() << std::endl;

	return 0;
}
