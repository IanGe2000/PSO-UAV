﻿# CMakeList.txt : Top-level CMake project file, do global configuration
# and include sub-projects here.
#
cmake_minimum_required (VERSION 3.8)

project ("PSO-UAV")

# Include sub-projects.
add_subdirectory ("CXX_simulation")
add_subdirectory ("PSO-UAV_Input")
