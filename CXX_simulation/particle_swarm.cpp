#include "particle_swarm.h"

RowVectorXd particleAdjust(RowVectorXd particle)
{
	//std::cout << "In: " << particle << "\n\n";
	Index i, k;
	double v = (diff(particle).array() >= 0).cast<double>().minCoeff(&k, &i);
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
	{
		//std::cout << "Out(Done): " << particle << "\n\n";
		return particle;
	}
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

Tensor<double, 3> velocityInit(RowVector2d velocity_range, int n, int N, int G_n)
{
	Tensor<double, 3> velocity(G_n, n, N);
	velocity.setRandom();
	double velocity_range_size = diff(velocity_range)(0);
	velocity = velocity * (velocity_range_size) + velocity_range.minCoeff();
	return velocity;
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
				objective(i, j, 0) = std::numeric_limits<double>::infinity();
			else
				objective(i, j, 0) = omega_d * F_d(trajectory_m) + omega_c * F_c(particle, i, solution, CT) + (1 - omega_c - omega_d) * F_a(trajectory_m, N_W);
		}
	}

	return objective;
}

int PSOMain(int n, int N, int G_n, Matrix2d startendpoints, double theta_Tmax, double theta_Cmax, Matrix3Xd threat_source, int maxgeneration, std::fstream* Log, std::fstream* Solution)
{
	// Step 1 //
	std::clock_t t0 = std::clock();
	bool success = true;
	bool timer_running = true;
	//int n = 10;		//number of particles in a subgroup
	//int N = 10;		//number of way-points
	//int G_n = 5;	//number of subgroups
	//Matrix2d startendpoints{ {0, 30},{0, 30} };
	Map<Matrix<double, 2, 1>> startpoint(startendpoints(placeholders::all, 0).data(), 2, 1);
	Map<Matrix<double, 2, 1>> endpoint(startendpoints(placeholders::all, 1).data(), 2, 1);
	RowVector2d position_range;
	position_range(0) = startendpoints(1, placeholders::all).minCoeff();
	position_range(1) = startendpoints(1, placeholders::all).maxCoeff();
	RowVectorXd xIntervals = VectorXd::LinSpaced(N, startendpoints(0, 0), startendpoints(0, 1));
	// swarm initialization
	Tensor<double, 3> swarm = swarmInit(position_range, n, N, G_n);
	// constraints
	//double theta_Tmax = 60;
	//double theta_Cmax = 45;
	//Matrix3Xd threat_source
	//{
	//	{6.2280, 17.781, 15.681, 6.5280, 22.581, 15.057, 21.036},
	//	{8.5230, 4.6080, 17.208, 13.629, 21.108, 11.835, 15.846},
	//	{2.2826, 1.9663, 2.8540, 2.0762, 1.9393, 2.4483, 2.4404}
	//};
	//int maxgeneration = 30;
	double P_c = 0.85;
	ArrayXd omega = ArrayXd::LinSpaced(maxgeneration, 0.7, 0.4);
	double phi_p = 0.2;    // cognitive coefficient
	double phi_g = 0.2;    // social coefficient
	RowVector2d velocity_range{ {-6.0, 6.0} };
	Tensor<double, 3> velocity = velocityInit(velocity_range, n, N, G_n);
	// objective function weights
	double omega_d = 0.3;
	double omega_c = 0.5;
	double CT = 4.5;
	int N_W = N;

	// Step 2 //
	int generation = 1;
	// Initialize the particle's best position to be the current swarm and the corresponding objective value to be inf
	Tensor<double, 3> P_pos = swarm;
	Tensor<double, 3> P_obj(G_n, n, 1);
	P_obj.setConstant(std::numeric_limits<double>::infinity());
	// Initialize the group best particle to be the 1st particle in its subgroup
	array<Index, 3> offsets_3_swarm = { 0,0,0 };
	array<Index, 3> extents_3_particle = { 1,1,N };
	array<Index, 3> extents_3_one_particle_of_each_group = { G_n,1,N };
	array<Index, 3> offsets_3_G_pos = { 0,0,0 };
	array<Index, 3> extents_3_group = { 1,n,N };
	array<Index, 3> extents_3_one_col_of_each_group = { G_n,n,1 };
	Tensor<double, 3> G_pos = swarm.slice(offsets_3_swarm, extents_3_one_particle_of_each_group);
	Tensor<double, 3> G_obj(G_n, 1, 1);
	G_obj.setConstant(std::numeric_limits<double>::infinity());
	// The solution is apparently the best particle from each group
	Map<MatrixXd> solution(G_pos.data(), G_n, N);

	// visualization
	*Log << "Initial generation:\n";
	*Log << "swarm:\n" << swarm << "\n\n";
	*Log << "solution:\n" << solution << "\n\n";
	*Log << "G_obj:\n" << G_obj << "\n\n";

	MatrixXd G_obj_log(G_n, maxgeneration);
	Map<MatrixXd> G_obj_temp(G_obj.data(), G_n, 1);
	G_obj_log(placeholders::all, 0) = G_obj_temp;

	while (generation <= maxgeneration)
	{
		*Log << "enter generation " << generation << ":================================================================================\n";

		// calculate the objective of this iteration of swarm
		Tensor<double, 3> P_objective = F(swarm, xIntervals, solution, threat_source, theta_Tmax, theta_Cmax, omega_d, omega_c, CT, N_W);
		// refresh P_posand P_obj
		Tensor<bool, 3> P_refresh = P_objective < P_obj;
		for (int i = 0; i < G_n; i++)
		{
			for (int j = 0; j < n; j++)
			{
				if (P_refresh(i, j, 0))
				{
					// particle's best is refreshed
					offsets_3_swarm = { i,j,0 };
					P_pos.slice(offsets_3_swarm, extents_3_particle) = swarm.slice(offsets_3_swarm, extents_3_particle);
					P_obj(i,j,0) = P_objective(i,j,0);
				}
			}
		}

		// Step 3 //
		// refresh G_pos and G_obj
		Tensor<int, 2> group_best_index = P_obj.argmin(1).cast<int>();
		for (int i = 0; i <	G_n; i++)
		{
			if (P_obj(i, group_best_index(i), 0) < G_obj(i, 0, 0))
			{
				// group best is refreshed
				offsets_3_swarm = { i, group_best_index(i), 0 };
				offsets_3_G_pos = { i, 0, 0 };
				G_pos.slice(offsets_3_G_pos, extents_3_particle) = swarm.slice(offsets_3_swarm, extents_3_particle);
				G_obj(i, 0, 0) = P_obj(i, group_best_index(i), 0);
			}
		}
		// regenerates the solution from G_pos
		new (&solution) Map<MatrixXd> (G_pos.data(), G_n, N);
		// check if there is any duplication of particles
		ArrayXi reinit_group_indexes = repeatedRow(solution);
		swarm = swarmInit(swarm, reinit_group_indexes);
		for (int i = 0; i < reinit_group_indexes.rows(); i++)
		{
			offsets_3_swarm = { reinit_group_indexes(i), 0, 0 };
			P_pos.slice(offsets_3_swarm, extents_3_group) = swarm.slice(offsets_3_swarm, extents_3_group);
			P_obj(reinit_group_indexes(i), 0, 0) = std::numeric_limits<double>::infinity();
			G_pos.slice(offsets_3_swarm, extents_3_particle) = swarm.slice(offsets_3_swarm, extents_3_particle);
			G_obj(reinit_group_indexes(i), 0, 0) = std::numeric_limits<double>::infinity();
		}

		new (&G_obj_temp) Map<MatrixXd>(G_obj.data(), G_n, 1);
		G_obj_log(placeholders::all, generation - 1) = G_obj_temp;

		// visualization
		*Log << "solution:\n" << solution << "\n\n";
		*Log << "G_obj:\n" << G_obj << "\n\n";
		*Log << "G_obj_log:\n" << G_obj_log << "\n\n";
		*Log << "NOW update velocity and position:------------------------------------------------------------------\n";

		// Step 4 //
		// update velocity and position of the swarm
		// modification: uses random C/D switching PSO with convergence ratio P_c for velocity update
		for (int i = 0; i < G_n; i++)
		{
			*Log << "Group " << i << ": ";
			if ((reinit_group_indexes == 2).any())
				continue;
			double xi = Array2d::Random().abs()(0);
			// force operator D when no feasable solution is found in this group
			if (G_obj(i, 0, 0) == std::numeric_limits<double>::infinity())
				xi = 2;
			if (xi <= P_c)
				*Log << "Operator C\n";
			else if (xi == 2)
				*Log << "Operator D (forced)\n";
			else
				*Log << "Operator D\n";
			for (int j = 0; j < n; j++)
			{
				for (int k = 0; k < N; k++)
				{
					double r1 = Array2d::Random().abs()(0);
					double r2 = Array2d::Random().abs()(0);
					if (xi <= P_c)
					{
						// Operator C
						velocity(i, j, k) = omega(generation - 1) * velocity(i, j, k) + phi_p * r1 * (P_pos(i, j, k) - swarm(i, j, k)) + (0.5 * r2 + 0.5) * (2 * omega(generation - 1) + 2 - phi_g * r2) * (G_pos(i, 0, k) - swarm(i, j, k));
					}
					else
					{
						// Operator D
						velocity(i, j, k) = 1.05 * velocity(i, j, k) + phi_p * r1 * (P_pos(i, j, k) - swarm(i, j, k)) + phi_g * r2 * (G_pos(i, 1, k) - swarm(i, j, k));
					}
					swarm(i, j, k) = swarm(i, j, k) + velocity(i, j, k);
					// submit to the start - end point and position range constraint
					if (k == 0)
						swarm(i, j, k) = startpoint(1);
					else if (k == N - 1)
						swarm(i, j, k) = endpoint(1);
					if (swarm(i, j, k) < position_range(0))
						swarm(i, j, k) = position_range(0);
					else if (swarm(i, j, k) > position_range(1))
						swarm(i, j, k) = position_range(1);
				}
			}
		}

		// Step 5 //
		// adjust the particles to satisfy y(i) <= y(i+1)
		for (int i = 0; i < G_n; i++)
		{
			//Log << "i = " << i << ":\n";
			for (int j = 0; j < n; j++)
			{
				//Log << "j = " << j << ":\n";
				offsets_3_swarm = { i,j,0 };
				Tensor<double, 3> swarm_temp(1, 1, N);
				swarm_temp = swarm.slice(offsets_3_swarm, extents_3_particle);
				//Log << "swarm_temp: " << swarm_temp << "\n";
				Map<RowVectorXd> particle_T2M(swarm_temp.data(), 1, N);
				//Log << "particle_T2M: " << particle_T2M << "\n";
				RowVectorXd result = particleAdjust(particle_T2M);
				//Log << "result: " << result << "\n";
				TensorMap<Tensor<double, 3>> particle_M2T(result.data(), 1, 1, N);
				//Log << "particle_M2T: " << particle_M2T << "\n";
				swarm.slice(offsets_3_swarm, extents_3_particle) = particle_M2T;
			}
		}

		// visualization
		*Log << "swarm:\n" << swarm << "\n\n";

		// timer
		if (timer_running)
		{
			for (int i = 0; i < G_n; i++)
			{
				if (G_obj(i, 0, 0) == std::numeric_limits<double>::infinity())
					success = success && false;
				else
					success = success && true;
			}
			if (success)
			{
				std::clock_t t1 = std::clock();
				double time = (t1 - t0) / (CLOCKS_PER_SEC / 1000);
				std::cout << time << "ms, at generation " << generation << "\n";
				timer_running = false;
			}
			else
				success = true;
		}
		generation++;
	}

	// visualization
	*Solution << G_n << "\n";
	*Solution << N << "\n";
	*Solution << xIntervals << "\n";
	*Solution << solution << "\n";

	return 0;
}

Matrix2Xd plotCircle(Matrix<double, 3, 1> threat_source)
{
	Array<double, 1, Dynamic> t = Array<double, 1, Dynamic>::LinSpaced(100, 0, 2 * acos(-1));
	Matrix2Xd circbmp(2, 100);
	circbmp(0, placeholders::all) = threat_source(2, 0) * cos(t) + threat_source(0, 0);
	circbmp(1, placeholders::all) = threat_source(2, 0) * sin(t) + threat_source(1, 0);
	return circbmp;
}

void plotThreat_source(Matrix3Xd threat_source, std::fstream* Threat_source_bmp)
{
	*Threat_source_bmp << threat_source.cols() << "\n";
	for (int i = 0; i < threat_source.cols(); i++)
	{
		*Threat_source_bmp << plotCircle(threat_source(placeholders::all, i)) << "\n";
	}
}
