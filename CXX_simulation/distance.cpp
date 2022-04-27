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

	return distance((MatrixXd)X(placeholders::all, seq(0, placeholders::last-1)), (MatrixXd)X(placeholders::all, seq(1, placeholders::last)));
}

MatrixXd distance(MatrixXd X, MatrixXd Y, distance_opt opt)
{
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
		{
			MatrixXd result(X.cols(), Y.cols());
			for (int i = 0; i < result.rows(); i++)
			{
				for (int j = 0; j < result.cols(); j++)
				{
					result(i, j) = (X(placeholders::all, i) - Y(placeholders::all, j)).array().square().colwise().sum().sqrt()(0);
				}
			}
			return result;
		}
		break;
	case recursive:
		if (X.rows() != 2 || Y.rows() != 2)
			errorHandler(only_point_to_line_distance_supported);
		else
		{
			MatrixXd result(Y.cols(), X.cols() - 1);
			Matrix2d num;
			Array<double, 1, 1> det;

			for (int i = 0; i < result.rows(); i++)
			{
				for (int j = 0; j < result.cols(); j++)
				{
					num << X(0, j) - X(0, j + 1), X(1, j) - X(1, j + 1), Y(0, i) - X(0, j + 1), Y(1, i) - X(1, j + 1);
					det << num.determinant();
					result(i, j) = det.abs()(0) / distance((MatrixXd)X(placeholders::all, j), (MatrixXd)X(placeholders::all, j + 1))(0);
				}
			}
			return result;
		}
		break;
	}
}

MatrixXd distance(MatrixXd X, MatrixXd A, MatrixXd B, distance_opt opt)
{
	switch (opt)
	{
	case parallel:
		if (X.rows() != A.rows() || X.rows() != B.rows() || X.cols() != A.cols() || X.cols() != B.cols())
			errorHandler(dimensions_mismatch);
		else if (X.rows() != 2)
			errorHandler(only_point_to_line_distance_supported);
		else
		{
			MatrixXd result(1, X.cols());
			Matrix2d num;
			Array<double, 1, 1> det;

			for (int i = 0; i < result.cols(); i++)
			{
				num << A(0, i) - B(0, i), A(1, i) - B(1, i), X(0, i) - B(0, i), X(1, i) - B(1, i);
				det << num.determinant();
				result(i) = det.abs()(0) / distance((MatrixXd)A(placeholders::all, i), (MatrixXd)B(placeholders::all, i))(0);
			}
			return result;
		}
		break;
	case mesh:
		if (A.rows() != B.rows() || A.cols() != B.cols())
			errorHandler(dimensions_mismatch);
		else if (X.rows() != 2)
			errorHandler(only_point_to_line_distance_supported);
		else
		{
			MatrixXd result(X.cols(), A.cols());
			Matrix2d num;
			Array<double, 1, 1> det;

			for (int i = 0; i < result.rows(); i++)
			{
				for (int j = 0; j < result.cols(); j++)
				{
					num << A(0, j) - B(0, j), A(1, j) - B(1, j), X(0, i) - B(0, j), X(1, i) - B(1, j);
					det << num.determinant();
					result(i) = det.abs()(0) / distance((MatrixXd)A(placeholders::all, j), (MatrixXd)B(placeholders::all, j))(0);
				}
			}
			return result;
		}
		break;
	default:
		errorHandler(unsupported_opt);
		break;
	}
}
