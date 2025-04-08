cmake_minimum_required(VERSION 3.1)

project(PQPIKS VERSION 1.0)

#Use C++11
set (CMAKE_CXX_STANDARD 11)

#Require Eigen 3
#find_package(Eigen 3 REQUIRED)

#Bring the qpOASES library into the project
set(PROJECT_LINK_LIBS libqpOASES.so)
link_directories(../qpOASES_master/bin)

#Bring the headers into the project (project headers and qpOASES headers)
include_directories(include
					../qpOASES_master/include)

#Generate the shared library from the sources
add_library(${PROJECT_NAME} SHARED
			src/PQPIKS.cpp
			src/PQPIKS_MDH.cpp
			src/KinematicModelMDH.cpp
			src/DynamicModel.cpp
			src/DynamicModelMDH.cpp
			src/Trajectory.cpp)

target_link_libraries(${PROJECT_NAME} ${PROJECT_LINK_LIBS})
target_compile_definitions(${PROJECT_NAME} PRIVATE __USE_LONG_INTEGERS__ __USE_LONG_FINTS__) # without it qpOASES is giving errors during compilation
