#include "distance.h"

std::string error_decoder[] = { "unsupported_opt", "dimensions_mismatch", "only_point_to_line_distance_supported"};

void errorHandler(distance_error error)
{
	std::cout << "distance: " << error_decoder[error] << std::endl;
	exit(1);
}

RowVectorXd distance(MatrixXd X, distance_opt opt)
{
	if (opt != recursive)
		errorHandler(unsupported_opt);

	return distance(X(placeholders::all, seq(0, placeholders::last-1)), X(placeholders::all, seq(1, placeholders::last)));
}

MatrixXd distance(MatrixXd X, MatrixXd Y, distance_opt opt)
{
	std::cout << X << "\n\n" << Y << std::endl;
	switch (opt)
	{
	case parallel:
		if (X.rows() != Y.rows() || X.cols() != Y.cols())
			errorHandler(dimensions_mismatch);
		else
			return (X - Y).array().square().colwise().sum().sqrt();
		break;
	case mesh:
		if (X.rows() != Y.rows())
			errorHandler(dimensions_mismatch);
		else
			// TODO mesh calculation
			return (X - Y).array().square().colwise().sum().sqrt();
		break;
	case recursive:
		if (X.rows() != 2 || Y.rows() != 2)
			errorHandler(only_point_to_line_distance_supported);
		else
			//TODO recursive calculation
			return (X - Y).array().square().colwise().sum().sqrt();
		break;
	}
}
