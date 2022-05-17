#include "CXX_simulation.h"
#include "utilities.h"
#include "altitude.h"

RowVectorXd particleAdjust(RowVectorXd);
Tensor<double, 3> swarmInit(RowVector2d, int, int, int);
Tensor<double, 3> swarmInit(Tensor<double, 3>, int);
Tensor<double, 3> swarmInit(Tensor<double, 3>, ArrayXi);
Tensor<double, 3> velocityInit(RowVector2d, int, int, int);
Matrix3Xd particle2Trajectory(RowVectorXd, RowVectorXd);
Tensor<double, 4> swarm2Trajectory(Tensor<double, 3>, RowVectorXd);
Matrix2Xd particle2Course(RowVectorXd, RowVectorXd);
double F_d(Matrix3Xd);
double F_d(RowVectorXd, RowVectorXd);
bool threatConflict(Matrix2Xd, Matrix3Xd);
double turningAngle(Matrix3d);
Array<bool, 1, Dynamic> maxTurningAngle(Matrix3Xd, double);
int F_c(RowVectorXd, int, MatrixXd, double);
double F_a(Matrix3Xd, int);
double climbingDivingAngle(Matrix2d);
Array<bool, 1, Dynamic> maxClimbingDivingAngle(Matrix2Xd, double);
Tensor<double, 3> F(Tensor<double, 3>, RowVectorXd, MatrixXd, Matrix3Xd, double, double, double, double, double, int);
int PSOMain(int, int, int, Matrix2d, double, double, Matrix3Xd, int, std::fstream*, std::fstream*);
Matrix2Xd plotCircle(Matrix<double, 3, 1>);
void plotThreat_source(Matrix3Xd, std::fstream*);
