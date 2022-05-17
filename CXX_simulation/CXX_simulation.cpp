#include "CXX_simulation.h"
#include "utilities.h"
#include "particle_swarm.h"
#include "altitude.h"

int main()
{
	// File Interface
	std::fstream Input_parameters("Input_parameters.txt", std::fstream::in);
	std::fstream Log("Log.txt", std::fstream::out | std::fstream::trunc);
	std::fstream Solution("Solution.txt", std::fstream::out | std::fstream::trunc);
	std::fstream Threat_source_bmp("Threat_source_bmp.txt", std::fstream::out | std::fstream::trunc);

	// Parameters declaration
	int G_n = 0;
	Matrix2d startendpoints;
	int number_of_threat_sources = 0;
	Matrix3Xd threat_source;
	int n = 0;
	int N = 0;
	double theta_Tmax = 0;
	double theta_Cmax = 0;
	int maxgeneration = 0;

	// Extract parameters from file
	std::string filecontent;
	std::stringstream ss;
	int switcher0 = 0;
	int switcher1 = 0;
	double temp = 0;
	while (getline(Input_parameters, filecontent))
	{
		switch (switcher0)
		{
		case 0:
			G_n = std::stoi(filecontent);
			switcher0++;
			break;
		case 1:
			ss << filecontent;
			for (int i = 0; i < 2; i++)
			{
				ss >> temp;
				startendpoints(i, 0) = temp;
			}
			switcher0++;
			break;
		case 2:
			new (&ss) std::stringstream;
			ss << filecontent;
			for (int i = 0; i < 2; i++)
			{
				ss >> temp;
				startendpoints(i, 1) = temp;
			}
			switcher0++;
			break;
		case 3:
			number_of_threat_sources = std::stoi(filecontent);
			new (&threat_source) Matrix3Xd(3, number_of_threat_sources);
			switcher0++;
			break;
		case 4:
			new (&ss) std::stringstream;
			ss << filecontent;
			for (int i = 0; i < number_of_threat_sources; i++)
			{
				ss >> threat_source(switcher1, i);
			}
			if (switcher1 < 2)
				switcher1++;
			else
			{
				switcher1 = 0;
				switcher0++;
			}
			break;
		case 5:
			n = std::stoi(filecontent);
			switcher0++;
			break;
		case 6:
			N = std::stoi(filecontent);
			switcher0++;
			break;
		case 7:
			new (&ss) std::stringstream;
			ss << filecontent;
			ss >> theta_Tmax >> theta_Cmax;
			switcher0++;
			break;
		case 8:
			maxgeneration = std::stoi(filecontent);
			switcher0++;
			break;
		default:
			if (filecontent.empty())
				break;
			else
			{
				std::cout << "Unexpected extra content detected.\n\n";
				break;
			}
		}
	}

	Input_parameters.close();

	PSOMain(n, N, G_n, startendpoints, theta_Tmax, theta_Cmax, threat_source, maxgeneration, &Log, &Solution);
	Log.close();
	Solution.close();

	plotThreat_source(threat_source, &Threat_source_bmp);
	Threat_source_bmp.close();
	return 0;
}
