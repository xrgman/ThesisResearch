cmake_minimum_required(VERSION 3.16.3)

project(demo)

# find_package(SimpleITK)
# add_executable ( sitk_example sitk_example.cxx )
set(CMAKE_CXX_STANDARD 20)
set(COMPILE_DEFINITIONS -O2 -Wall)

#Linking to cmake file that locates PortAudio
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}")

# Include the CMakeLists.txt in the "lib" folder
add_subdirectory(lib)

# Find the PortAudio package
find_package(PortAudio REQUIRED) 

#Creating executable
add_executable( demo src/main.cpp )

# Include PortAudio headers
target_include_directories(demo PRIVATE ${PORTAUDIO_INCLUDE_DIRS})

# Linking libraries
target_link_libraries(demo PRIVATE ${PORTAUDIO_LIBRARIES})
target_link_libraries(demo PRIVATE WavReader)
target_link_libraries(demo PRIVATE AudioHelper)
