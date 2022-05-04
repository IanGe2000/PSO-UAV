#pragma once

#include <string>
#include "CXX_simulation.h"

enum distance_opt
{
	parallel, mesh, recursive
};

enum utilities_error
{
	unsupported_opt, dimensions_mismatch, only_point_to_line_distance_supported, K_must_be_positive, DIM_must_be_1_or_2
};

void errorHandler(utilities_error);

RowVectorXd distance(MatrixXd, distance_opt = recursive);
MatrixXd distance(MatrixXd, MatrixXd, distance_opt = parallel);
MatrixXd distance(MatrixXd, MatrixXd, MatrixXd, distance_opt = parallel);

MatrixXd diff(MatrixXd);
MatrixXd diff(MatrixXd, int);
MatrixXd diff(MatrixXd, int, int);

