# CMakeList.txt : CMake project for CXX_simulation, include source and define
# project specific logic here.
#
cmake_minimum_required (VERSION 3.8)

# Set the include path
include_directories("../..")

# Add source to this project's executable.
add_executable (CXX_simulation "CXX_simulation.cpp" "CXX_simulation.h" "distance.cpp" "distance.h")

# TODO: Add tests and install targets if needed.
