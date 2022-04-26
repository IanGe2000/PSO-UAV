#pragma once

#include <iostream>

//#define USER_EIGEN

#ifdef USER_EIGEN
#include <Eigen/Dense>
#else
#include "Eigen/Dense"
#endif