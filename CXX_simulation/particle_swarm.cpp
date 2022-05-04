#include "particle_swarm.h"

RowVectorXd particleAdjust(RowVectorXd particle)
{
	std::cout << "In: " << particle << "\n\n";
	Index i, k;
	double v = (diff(particle).array() > 0).cast<double>().minCoeff(&k, &i);
	std::cout << "i = " << i << "\n\n";
	if (v == 0)
	{
		RowVectorXd::Index j;
		double w = ((particle(placeholders::all, seq(i, placeholders::last)).array() - particle(i)) > 0).cast<int>().maxCoeff(&k, &j);
		std::cout << "j = " << j << "\n\n";
		if (w == 1)
		{
			double ran = (Array<double, 1, 1>::Random()(0) / 2 + 0.5);
			std::cout << "ran = " << ran << "\n\n";
			particle(i + 1) = particle(i) + (particle(i + j) - particle(i)) / (j - 1) * ran;
			std::cout << "Out:" << particle << "\n\n";
			return particleAdjust(particle);
		}
		else
		{
			particle(seq(i,placeholders::last - 1)).array() = particle(i);
			std::cout << "Out:" << particle << "\n\n";
			return particleAdjust(particle);
		}
	}
	else
		return particle;
}
