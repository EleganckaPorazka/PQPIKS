cmake_minimum_required(VERSION 3.1)

project(example VERSION 1.0)

#Use C++11
set (CMAKE_CXX_STANDARD 11)

#Bring the PQPIKS and qpOASES libraries into the project
#(yes, qpOASES was already added to PQPIKS, but it seems that it has to be added again to be seen during compilation)
set(PROJECT_LINK_LIBS libPQPIKS.so
						libqpOASES.so)
link_directories(../build
				../../qpOASES_master/bin)

#Bring the headers into the project (InvKinMPC and qpOASES headers)
include_directories(../include
					../../qpOASES_master/include)

#Generate the executable from the source file
add_executable(example example.cpp)

#Link the libraries to the executable
target_link_libraries(example ${PROJECT_LINK_LIBS})
target_compile_definitions(${PROJECT_NAME} PRIVATE __USE_LONG_INTEGERS__ __USE_LONG_FINTS__) # without it qpOASES is giving errors during compilation
