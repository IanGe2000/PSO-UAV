#pragma once

#include <fstream>
#include <iostream>
#include <ctime>

#define USER_EIGEN

#ifdef USER_EIGEN
#include <Eigen/Dense>
#include <Eigen/unsupported/CXX11/Tensor>
#else
#include "Eigen/Dense"
#include "Eigen/unsupported/CXX11/Tensor"
#endif

using namespace Eigen;

