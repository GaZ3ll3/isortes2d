cmake_minimum_required(VERSION 2.8)
project(isortes)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS_RELEASE} -Ofast -DRUN_OMP -march=native -std=c++11 -fopenmp ")

link_directories("${OPENBLAS_PATH}/lib")
set(SOURCE_FILES main.cpp bbfmm.cpp blas_wrapper.cpp gmres.cpp geometry.cpp modes.cpp polynomials.cpp)
add_executable(isortes ${SOURCE_FILES})
target_link_libraries(isortes m blas triangle quadmath)

