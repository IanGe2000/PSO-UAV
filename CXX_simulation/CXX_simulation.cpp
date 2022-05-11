#include "CXX_simulation.h"
#include "utilities.h"
#include "particle_swarm.h"
#include "altitude.h"

int main()
{
	Matrix<double, 2, 2> a;
	a << 1,2.3, 2, 3;
	Array<bool, 2, 2> b;
	b = a.array() < 3.01;
	std::cout << !b.all();
	return 0;
}
