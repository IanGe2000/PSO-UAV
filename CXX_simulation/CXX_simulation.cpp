#include "CXX_simulation.h"
#include "distance.h"

using namespace Eigen;

int main()
{
	MatrixXd X(2,3);
	X << 0, 2, 0, 0, 0, 2;
	std::cout << X << "\n\n";
	std::cout << X(placeholders::all, seq(0, placeholders::last - 1)) << "\n\n" << X(placeholders::all, seq(1, placeholders::last)) << "\n\n";
	RowVectorXd r1 = distance(X);
	std::cout << r1 << std::endl;
	return 0;
}
