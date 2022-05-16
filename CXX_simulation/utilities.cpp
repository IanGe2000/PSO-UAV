#include "utilities.h"

std::string error_decoder[] = { "unsupported_opt", "dimensions_mismatch", "only_point_to_line_distance_supported", "K_must_be_positive", "DIM_must_be_1_or_2", "index_out_of_bound"};

void errorHandler(utilities_error error)
{
	std::cout << "utilities/errorHandler: " << error_decoder[error] << std::endl;
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
	//std::cout << "X: " << X << "\n\nY: " << Y << "\n\n";
	switch (opt)
	{
	case parallel:
		//std::cout << "parallel\n";
		if (X.rows() != Y.rows() || X.cols() != Y.cols())
			errorHandler(dimensions_mismatch);
		else {
			//std::cout << (X - Y).array().square().colwise().sum().sqrt() << "\n\n";
			return (X - Y).array().square().colwise().sum().sqrt();
		}
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
		//std::cout << "recursive\n";
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
					//std::cout << "num: " << num << "\n\n";
					det << num.determinant();
					//std::cout << "det: " << det << "\n\n";
					result(i, j) = det.abs()(0) / distance((MatrixXd)X(placeholders::all, j), (MatrixXd)X(placeholders::all, j + 1))(0);
				}
			}
			//std::cout << result << "\n\n";
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

MatrixXd diff(MatrixXd X)
{
	return diff(X, 1);
}

MatrixXd diff(MatrixXd X, int K)
{
	if (X.rows() > K)
		return diff(X, K, 1);
	else if (X.cols() > K - X.rows() + 1)
		return diff(X, K - X.rows() + 1, 2);
	else
	{
		Matrix<double, 0, 0> N;
		return N;
	}
}

MatrixXd diff(MatrixXd X, int K, int DIM)
{
	if (K == 0)
		return X;
	else if (K >= 1)
	{
		if (DIM == 1 && X.rows() > 1)
		{
			MatrixXd D(X.rows() - 1, X.cols());
			D = X(seq(1, placeholders::last), placeholders::all) - X(seq(0, placeholders::last - 1), placeholders::all);
			return diff(D, K - 1, DIM);
		}
		else if (DIM == 2 && X.cols() > 1)
		{
			MatrixXd D(X.rows(), X.cols() - 1);
			D = X(placeholders::all, seq(1, placeholders::last)) - X(placeholders::all, seq(0, placeholders::last - 1));
			return diff(D, K - 1, DIM);
		}
		else if (DIM != 1 && DIM != 2)
			errorHandler(DIM_must_be_1_or_2);
		else
		{
			Matrix<double, 0, 0> N;
			return N;
		}
	}
	else
		errorHandler(K_must_be_positive);
}

ArrayXi repeatedRow(MatrixXd solution)
{
	Matrix<bool, 1, Dynamic> row_status(1, solution.rows());
	Matrix<int, 1, Dynamic> repeated_list(1, solution.rows());
	row_status.setConstant(false);
	repeated_list.setConstant(-1);
	int list_length = 0;
	//std::cout << row_status << "\n";
	for (int i = 0; i < row_status.cols(); i++)
	{
		if (row_status(i))
			continue;
		//std::cout << "i = " << i << ":\n";
		for (int j = i + 1; j < row_status.cols(); j++)
		{
			if (row_status(j))
				continue;
			//std::cout << "j = " << j << ": ";
			if (solution(i, placeholders::all) == solution(j, placeholders::all))
			{
				//std::cout << " repeated, row_status = ";
				repeated_list(list_length++) = j;
				row_status(i) = true;
				row_status(j) = true;
				//std::cout << row_status;
			}
			//std::cout << "\n";
		}
	}
	Map<ArrayXi> repeated_list_sliced(repeated_list.data(), list_length);
	//std::cout << repeated_list << "\n\n";
	return repeated_list_sliced;
}
