cmake_minimum_required(VERSION 3.16.3)

project(demo)

# find_package(SimpleITK)
# add_executable ( sitk_example sitk_example.cxx )
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(COMPILE_DEFINITIONS -O2 -Wall)

#Linking to cmake file that locates PortAudio
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}")

# Include the CMakeLists.txt in the "lib" folder
add_subdirectory(lib)

# Find the PortAudio package
find_package(PortAudio REQUIRED) 

# Find Boost (required for gnuplot-iostream)
find_package(Boost COMPONENTS iostreams REQUIRED)
find_package(Threads REQUIRED)

#Creating executable
add_executable( demo src/main.cpp )

# Include PortAudio headers
target_include_directories(demo PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src ${PORTAUDIO_INCLUDE_DIRS} ${Boost_INCLUDE_DIRS} kissfft)

# Linking libraries
target_link_libraries(demo PRIVATE ${PORTAUDIO_LIBRARIES} ${Boost_LIBRARIES} WavHelper AudioHelper kissfft FFTWrapper Util) 

# Copy gnuplot-iostream.h to the build directory
#configure_file(src/gnuplot-iostream.h ${CMAKE_CURRENT_BINARY_DIR}/gnuplot-iostream.h COPYONLY)
