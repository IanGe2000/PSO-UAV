#include "CXX_simulation.h"
#include "utilities.h"
#include "particle_swarm.h"
#include "altitude.h"

int main()
{
	Tensor<double, 3> swarm(2, 2, 5);
	swarm.setValues({ {{0,0,0,2,3},{0,0,0,2,3}},{{0,0,0,2,3},{0,0,0,2,3}} });
	std::cout << swarm << "\n\n";

	RowVectorXd xInt(1, 5);
	xInt << 0, 1, 2, 3, 4;
	std::cout << xInt << "\n\n";

	Tensor<double, 4> trajectory = swarm2Trajectory(swarm, xInt);
	std::cout << trajectory << "\n\n";
	return 0;
}
