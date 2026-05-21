cmake_minimum_required(VERSION 3.10)
project(p4_main_sor)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Указываем путь к папке shared через абсолютный путь
# Выходим из p4_main_sor, выходим из core, заходим в shared
include_directories("${CMAKE_SOURCE_DIR}/../../shared")

add_executable(main main.cpp solver.cpp)