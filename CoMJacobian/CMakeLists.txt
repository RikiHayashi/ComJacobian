cmake_minimum_required(VERSION 2.8)

project(CoMJacobian)
set(CMAKE_BUILD_TYPE Debug)
find_package(Eigen3 REQUIRED)
find_package(Boost COMPONENTS system thread filesystem)
add_definitions(-lstdc++ -std=c++11)

include_directories(
	${BOOST_INCLUDE_DIRS}
	${Boost_INCLUDE_DIR}
	${EIGEN3_INCLUDE_DIR}
)

link_directories(
	${Boost_LIBRARY_DIRS}
	${Boost_LIBRARY_DIR}
)

add_executable(CoMJacobian Kinematics.cpp Link.cpp testCoMJacobian.cpp Kinematics.h Link.h ComJacobian.cpp ComJacobian.h Jacobian.cpp )

target_link_libraries(CoMJacobian ${Boost_LIBRARIES})
