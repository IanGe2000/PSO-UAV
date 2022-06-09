# CMakeList.txt : CMake project for PSO-UAV_Input, include source and define
# project specific logic here.
#
cmake_minimum_required (VERSION 3.8)

set(CMAKE_INCLUDE_CURRENT_DIR ON)

set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Qt5 COMPONENTS Widgets Charts)

set(project_sources
    main.cpp
    mainwindow.cpp
    threatsource.cpp
    solutionchart.cpp
    )
set(project_headers
    mainwindow.h
    threatsource.h
    solutionchart.h
    )

add_executable(PSO-UAV_Input ${project_sources} ${project_headers})

target_link_libraries(PSO-UAV_Input
    PUBLIC
    Qt5::Widgets
    Qt5::Charts)
