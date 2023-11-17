cmake_minimum_required(VERSION 3.16.3)

project(demo)

# find_package(SimpleITK)
# add_executable ( sitk_example sitk_example.cxx )
set(CMAKE_CXX_STANDARD 20)
set(COMPILE_DEFINITIONS -O2 -Wall)

add_executable( demo src/main.cpp )