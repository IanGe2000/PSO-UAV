#include "CXX_simulation.h"
#include "utilities.h"
#include "particle_swarm.h"

int main()
{
	RowVectorXd M(10), N;
	M << 1, 2, 3, 4, 6, 5, 3, 7, 8, 10;
	//N = diff(M,1,2);
	N = particleAdjust(M);
	std::cout << M <<"\n\n" << N;
	return 0;
}
