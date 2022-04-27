#pragma once

#include <iostream>

//#define USER_EIGEN

#ifdef USER_EIGEN
#include <Eigen/Dense>
#include <Eigen/unsupported/CXX11/Tensor>
#else
#include "Eigen/Dense"
#include "Eigen/unsupported/CXX11/Tensor"
#endif

using namespace Eigen;
