#include "CXX_simulation.h"
#include "distance.h"

using namespace Eigen;

int main()
{
	MatrixXd X(2, 2);
	X << 0, 0, 0, 2;
	MatrixXd Y(4,1);
	Y << 6, 3,0,0;
	std::cout << X << "\n\n" << Y(2) << "\n\n";
	//Matrix2d Z;
	//Z << Y(0, 0), X(0, 1), Y(1, 2), X(1, 2);
	//Array<double,1,1> Zd;
	//Zd << Z.determinant();
	//std::cout << Zd << std::endl;
	//std::cout << Zd.abs() << std::endl;
	std::cout << X.size() << "\n\n" << Y.size() << "\n\n";
	if (X.size() == Y.size())
		std::cout << "same" << "\n\n";
	else
		std::cout << "notsame" << "\n\n";
	//RowVectorXd r1 = distance(X);
	//std::cout << r1 << std::endl;
	return 0;
}
