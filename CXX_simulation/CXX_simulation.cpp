// CXX_simulation.cpp : Defines the entry point for the application.
//

#include "CXX_simulation.h"

//using namespace std;

using namespace Eigen;

//using Eigen::ArrayXXd;
//using Eigen::MatrixXd;
//using Eigen::seq;
//using Eigen::placeholders::all;

int main()
{
	MatrixXd m = MatrixXd::Random(4, 4);
	//MatrixXd n(3, 3);
	//MatrixXd r1(4, 4);
	//MatrixXd r2(3, 3);
	//m << 1, 2, 3, 4, 5, 6, 7, 8, 9;
	//n << 1, 2, 3, 4, 5, 6, 7, 8, 9;
	//r1 = m(seq(2, placeholders::last));
	//r2 = r1 * r1;
	std::cout << m;
}