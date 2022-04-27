#pragma once

#include <string>
#include "CXX_simulation.h"

using namespace Eigen;

enum distance_opt
{
	parallel, mesh, recursive
};

enum distance_error
{
	unsupported_opt, dimensions_mismatch, only_point_to_line_distance_supported
};

void errorHandler(distance_error);
RowVectorXd distance(MatrixXd, distance_opt = recursive);
MatrixXd distance(MatrixXd, MatrixXd, distance_opt = parallel);
MatrixXd distance(MatrixXd, MatrixXd, MatrixXd, distance_opt = parallel);
