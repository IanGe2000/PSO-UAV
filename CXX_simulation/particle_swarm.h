#include "CXX_simulation.h"
#include "utilities.h"
#include "altitude.h"

RowVectorXd particleAdjust(RowVectorXd);
Tensor<double, 3> swarmInit(RowVector2d, int, int, int);
Tensor<double, 3> swarmInit(Tensor<double, 3>, int);
Tensor<double, 3> swarmInit(Tensor<double, 3>, ArrayXi);
Matrix3Xd particle2Trajectory(RowVectorXd, RowVectorXd);
Tensor<double, 4> swarm2Trajectory(Tensor<double, 3>, RowVectorXd);