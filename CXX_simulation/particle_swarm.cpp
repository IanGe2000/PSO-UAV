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
			std::cout << "particle: " << particle << "\n";
			Matrix3Xd trajectory_temp_m = particle2Trajectory(particle, xIntervals);
			TensorMap<Tensor<double, 4>> trajectory_temp_t(trajectory_temp_m.data(), 1, 1, trajectory_temp_m.rows(), trajectory_temp_m.cols());
			trajectory.slice(traj_offsets, traj_extents) = trajectory_temp_t;
		}
	}

	return trajectory;
}
