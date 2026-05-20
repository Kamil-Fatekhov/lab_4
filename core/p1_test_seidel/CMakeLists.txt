cmake_minimum_required(VERSION 3.10)
project(DirichletPoissonTBB)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Ищем Intel TBB
find_package(TBB REQUIRED)

add_executable(solver_app main.cpp solver.cpp)
target_link_libraries(solver_app PRIVATE TBB::tbb)