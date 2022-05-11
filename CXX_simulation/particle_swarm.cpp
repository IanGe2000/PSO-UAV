#include "particle_swarm.h"

RowVectorXd particleAdjust(RowVectorXd particle)
{
	//std::cout << "In: " << particle << "\n\n";
	Index i, k;
	double v = (diff(particle).array() > 0).cast<double>().minCoeff(&k, &i);
	//std::cout << "i = " << i << "\n\n";
	if (v == 0)
	{
		RowVectorXd::Index j;
		double w = ((particle(placeholders::all, seq(i, placeholders::last)).array() - particle(i)) > 0).cast<int>().maxCoeff(&k, &j);
		//std::cout << "j = " << j << "\n\n";
		if (w == 1)
		{
			double ran = (Array<double, 1, 1>::Random()(0) / 2 + 0.5);
			//std::cout << "ran = " << ran << "\n\n";
			particle(i + 1) = particle(i) + (particle(i + j) - particle(i)) / (j - 1) * ran;
			//std::cout << "Out:" << particle << "\n\n";
			return particleAdjust(particle);
		}
		else
		{
			particle(seq(i,placeholders::last - 1)).array() = particle(i);
			//std::cout << "Out:" << particle << "\n\n";
			return particleAdjust(particle);
		}
	}
	else
		return particle;
}

Tensor<double, 3> swarmInit(RowVector2d position_range, int n, int N, int G_n)
{
	Tensor<double, 3> swarm(G_n, n, N);
	swarm.setRandom();
	//std::cout << swarm << "\n\n";
	Tensor<double, 3> a(G_n, n, N);
	double position_range_size = diff(position_range)(0);
	swarm = (swarm * 0.85 - a.setConstant(0.425)) * position_range_size;
	//std::cout << swarm << "\n\n";
	for (int i = 0; i < N; i++)
	{
		double linear = position_range_size / (N - 1) * i + position_range(0);
		for (int j = 0; j < G_n; j++)
		{
			for (int k = 0; k < n; k++)
			{
				swarm(j, k, i) = swarm(j, k, i) + linear;
				if (swarm(j, k, i) < position_range(0) || i == 0)
					swarm(j, k, i) = position_range(0);
				else if (swarm(j, k, i) > position_range(1) || i == N - 1)
					swarm(j, k, i) = position_range(1);
			}
		}
	}
	return swarm;
}

Tensor<double, 3> swarmInit(Tensor<double, 3> swarm, int reinit_group_index)
{
	int G_n = swarm.dimension(0);
	if (reinit_group_index >= G_n)
		errorHandler(index_out_of_bound);
	int n = swarm.dimension(1);
	int N = swarm.dimension(2);
	RowVector2d position_range = { swarm(0,0,0), swarm(0,0,N-1) };

	array<Index, 3> offset = { reinit_group_index,0,0 };
	array<Index, 3> extent = { 1,n,N };

	swarm.slice(offset, extent) = swarmInit(position_range, n, N, 1);
	return swarm;
}

Tensor<double, 3> swarmInit(Tensor<double, 3> swarm, ArrayXi reinit_group_indexes)
{
	for (int i = 0; i < reinit_group_indexes.rows(); i++)
		swarm = swarmInit(swarm, reinit_group_indexes(i));
	return swarm;
}

Matrix3Xd particle2Trajectory(RowVectorXd particle, RowVectorXd xIntervals)
{
	Matrix2Xd course(2, particle.cols());
	course << xIntervals, particle;
	Matrix3Xd trajectory(3, particle.cols());
	trajectory << course, altitude(course);
	return trajectory;
}

Tensor<double, 4> swarm2Trajectory(Tensor<double, 3> swarm, RowVectorXd xIntervals)
{
	int G_n = swarm.dimension(0);
	int n = swarm.dimension(1);
	int N = swarm.dimension(2);

	Tensor<double, 4> trajectory(G_n, n, 3, N);
	array<Index, 4> traj_offsets = { 0,0,0,0 };
	array<Index, 4> traj_extents = { 1,1,3,N };
	Tensor<double, 3> swarm_temp(1, 1, N);
	array<Index, 3> swar_offsets = { 0,0,0 };
	array<Index, 3> swar_extents = { 1,1,N };

	for (int i = 0; i < G_n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			traj_offsets = { i,j,0,0 };
			swar_offsets = { i,j,0 };
			swarm_temp = swarm.slice(swar_offsets, swar_extents);
			Map<RowVectorXd> particle(swarm_temp.data(), 1, N);
			//std::cout << "particle: " << particle << "\n";
			Matrix3Xd trajectory_temp_m = particle2Trajectory(particle, xIntervals);
			TensorMap<Tensor<double, 4>> trajectory_temp_t(trajectory_temp_m.data(), 1, 1, trajectory_temp_m.rows(), trajectory_temp_m.cols());
			trajectory.slice(traj_offsets, traj_extents) = trajectory_temp_t;
		}
	}

	return trajectory;
}

Matrix2Xd particle2Course(RowVectorXd particle, RowVectorXd xIntervals)
{
	Matrix2Xd course(2, particle.cols());
	course << xIntervals, particle;
	return course;
}

double F_d(Matrix3Xd trajectory)
{
	double den = distance(trajectory(placeholders::all, seq(0, 0)), trajectory(placeholders::all, placeholders::last))(0);
	double num = distance(trajectory).sum();
	return num / den;
}

double F_d(RowVectorXd particle, RowVectorXd xIntervals)
{
	Matrix3Xd trajectory = particle2Trajectory(particle, xIntervals);
	return F_d(trajectory);
}

bool threatConflict(Matrix2Xd course, Matrix3Xd threat_source)
{
	int N = course.cols();
	int T = threat_source.cols();

	MatrixXd d1 = distance(course, threat_source(seq(0, 1), placeholders::all), recursive);
	MatrixXd d2 = distance(course, threat_source(seq(0, 1), placeholders::all), mesh);
	MatrixXd d3 = distance(course);
	//std::cout << d1 << "\n\n" << d2 << "\n\n" << d3 << "\n\n";

	Array<bool, Dynamic, Dynamic> booltable0(N, T);
	Array<bool, Dynamic, Dynamic> booltable1(T, N - 1);
	Array<bool, Dynamic, Dynamic> booltable2(T, N - 1);
	Array<bool, Dynamic, Dynamic> booltable3(T, N - 1);
	Array<bool, Dynamic, Dynamic> booltable(T, N - 1);

	for (int i = 0; i < N; i++)
		booltable0(i/**/, placeholders::all) = d2(i/**/, placeholders::all).array() < threat_source(2, placeholders::all).array();
	for (int i = 0; i < T; i++)
		booltable1(i/**/, placeholders::all) = (booltable0(seq(0, placeholders::last - 1), i/**/) || booltable0(seq(1, placeholders::last), i/**/)).transpose();
	for (int i = 0; i < N - 1; i++)
		booltable2(placeholders::all, i/**/) = d1(placeholders::all, i/**/).array() > threat_source(2, placeholders::all).transpose().array();
	for (int i = 0; i < T; i++)
		booltable3(i/**/, placeholders::all) = d3.array().square() > diff(d2.array().square()).array().abs().transpose()(i/**/, placeholders::all);

	booltable = booltable1 || !booltable1 && !booltable2 && booltable3;
	//std::cout << booltable0 << "\n\n" << booltable1 << "\n\n" << booltable2 << "\n\n" << booltable3 << "\n\n" << booltable << "\n\n";
	return booltable.maxCoeff();
}

double turningAngle(Matrix3d segment)
{
	Matrix<double, 3, 2> dif = diff(segment,1,2);
	double theta_rad = acos((dif(0) * dif(3) + dif(1) * dif(4)) / sqrt((pow(dif(0), 2) + pow(dif(1), 2)) * (pow(dif(3), 2) + pow(dif(4), 2))));
	double theta_deg = theta_rad / acos(-1) * 180;
	return theta_deg;
}

Array<bool, 1, Dynamic> maxTurningAngle(Matrix3Xd trajectory, double theta_Tmax)
{
	int N = trajectory.cols();
	Array<bool, 1, Dynamic> boolarray(1, N - 2);
	for (int i = 0; i < boolarray.cols(); i++)
		boolarray(i) = turningAngle(trajectory(placeholders::all, seq(i, i + 2))) < theta_Tmax;
	return boolarray;
}

int F_c(RowVectorXd particle, int groupindex, MatrixXd solution, double CT)
{
	MatrixXd D_M(solution.rows() - 1, solution.cols());
	int k = 0;
	for (int i = 0; i < solution.rows(); i++)
	{
		if (i == groupindex)
			continue;
		else
			D_M(k++, placeholders::all) = (particle - solution(i, placeholders::all)).array().abs();
	}
	return (D_M.colwise().minCoeff().array() < CT).cast<int>().sum();
}

double F_a(Matrix3Xd trajectory, int N_W)
{
	return trajectory(2, placeholders::all).sum() / trajectory.cols();
}

double climbingDivingAngle(Matrix2d slice)
{
	double theta_rad = atan(diff(altitude(slice)).array().abs()(0) / sqrt(pow(diff(slice, 2).array(), 2).sum()));
	double theta_deg = theta_rad / acos(-1) * 180;
	return theta_deg;
}

Array<bool, 1, Dynamic> maxClimbingDivingAngle(Matrix2Xd course, double theta_Cmax)
{
	int N = course.cols();
	Array<bool, 1, Dynamic> boolarray(1, N - 1);
	for (int i = 0; i < boolarray.cols(); i++)
		boolarray(i) = climbingDivingAngle(course(placeholders::all, seq(i, i + 1))) < theta_Cmax;
	return boolarray;
}

Tensor<double, 3> F(Tensor<double, 3> swarm, RowVectorXd xIntervals, MatrixXd solution, Matrix3Xd threat_source, double theta_Tmax, double theta_Cmax, double omega_d, double omega_c, double CT, int N_W)
{
	int G_n = swarm.dimension(0);
	int n = swarm.dimension(1);
	int N = swarm.dimension(2);
	Tensor<double, 3> objective(G_n, n, 1);
	Tensor<double, 4> trajectory = swarm2Trajectory(swarm, xIntervals);
	Tensor<double, 4> trajectory_temp(1, 1, 3, N);
	array<Index, 4> traj_offsets = { 0,0,0,0 };
	array<Index, 4> traj_extents = { 1,1,3,N };
	Tensor<double, 3> swarm_temp(1, 1, N);
	array<Index, 3> swar_offsets = { 0,0,0 };
	array<Index, 3> swar_extents = { 1,1,N };

	for (int i = 0; i < G_n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			swar_offsets = { i,j,0 };
			swarm_temp = swarm.slice(swar_offsets, swar_extents);
			Map<RowVectorXd> particle(swarm_temp.data(), 1, N);
			traj_offsets = { i,j,0,0 };
			trajectory_temp = trajectory.slice(traj_offsets, traj_extents);
			Map<Matrix3Xd> trajectory_m(trajectory_temp.data(), 3, N);
			if (threatConflict(particle2Course(particle, xIntervals), threat_source) || !maxTurningAngle(trajectory_m, theta_Tmax).all() || !maxClimbingDivingAngle(trajectory_m(seq(0,1), placeholders::all), theta_Cmax).all())
				objective(i, j, 0) = Infinity;
			else
				objective(i, j, 0) = omega_d * F_d(trajectory_m) + omega_c * F_c(particle, i, solution, CT) + (1 - omega_c - omega_d) * F_a(trajectory_m, N_W);
		}
	}
}
