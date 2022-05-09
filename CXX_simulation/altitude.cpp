#include "altitude.h"

RowVectorXd altitude(MatrixXd way_point)
{
	static double p1 = 0.1501;
	static double p2 = 0.2311;
	static double p3 = 0.7068;
	static double p4 = 0.4860;
	static double p5 = 0.6913;
	static double p6 = 0.2621;
	static double p7 = 0.4565;
	static Array<double, 1, Dynamic> x;
	x = way_point(seq(0, 0), placeholders::all).array();
	static Array<double, 1, Dynamic> y;
	y = way_point(seq(1, 1), placeholders::all).array();
	RowVectorXd altitude = sin(y + p1) + p2*sin(x) + p3*cos(p4*sqrt(x.square() + y.square())) + p5*cos(y) + p6*sin(p6*sqrt(x.square() + y.square())) + p7*cos(y);
	return altitude;
}