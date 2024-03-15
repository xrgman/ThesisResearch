cmake_minimum_required(VERSION 3.16.3)

project(demo)

# find_package(SimpleITK)
# add_executable ( sitk_example sitk_example.cxx )
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(COMPILE_DEFINITIONS -O2 -Wall)

#Linking to cmake file that locates PortAudio
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}")

#Find json library:
#include_directories(${CMAKE_CURRENT_SOURCE_DIR}/lib/json.hpp)
# find_package(nlohmann_json 3.2.0 REQUIRED)
# add_subdirectory(lib/nlohmann_json)
include(FetchContent)

#FetchContent_Declare(json URL https://github.com/nlohmann/json/releases/download/v3.11.3/json.tar.xz)
FetchContent_Declare(
  json
  GIT_REPOSITORY https://github.com/nlohmann/json.git
  GIT_TAG        v3.11.3
)
FetchContent_MakeAvailable(json) 

# Include the CMakeLists.txt in the "lib" folder
add_subdirectory(lib) 

# Find the PortAudio package
# find_package(PortAudio REQUIRED) 
# link_directories(${PORTAUDIO_LIBRARY_DIRS})

FetchContent_Declare(
  portaudio
  GIT_REPOSITORY https://github.com/PortAudio/portaudio.git
  GIT_TAG        master  # Or use a specific tag or commit hash
)

FetchContent_MakeAvailable(portaudio)

#Find SDL2
FetchContent_Declare(
  SDL2
  GIT_REPOSITORY https://github.com/libsdl-org/SDL.git
  GIT_TAG        main  # Or use a specific tag or commit hash
)

FetchContent_MakeAvailable(SDL2)

# Fetch SDL2_ttf
FetchContent_Declare(
  SDL2_ttf
  GIT_REPOSITORY https://github.com/libsdl-org/SDL_ttf.git
  GIT_TAG        main  # Or use a specific tag or commit hash
)
FetchContent_MakeAvailable(SDL2_ttf)
# find_package(SDL2 REQUIRED)   
# find_package(SDL2_ttf REQUIRED)
# include_directories(${SDL2_INCLUDE_DIRS} ${SDL2_TTF_INCLUDE_DIRS})

# Find Boost (required for gnuplot-iostream)
find_package(Boost COMPONENTS iostreams REQUIRED)
find_package(Threads REQUIRED)

# Link pthread library
find_package(Threads REQUIRED)

# Link spdlog
find_package(spdlog REQUIRED) 

#Creating executable
add_executable(demo src/main.cpp)  

# Include PortAudio headers
target_include_directories(demo PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src ${CMAKE_CURRENT_SOURCE_DIR}/lib/Config ${Boost_INCLUDE_DIRS} kissfft)
# target_include_directories(demo PRIVATE ${PORTAUDIO_INCLUDE_DIRS}) 
#target_include_directories(demo PRIVATE ${SDL2_INCLUDE_DIRS} ${SDL2_TTF_INCLUDE_DIRS})

# Linking libraries
target_link_libraries(demo PRIVATE ${Boost_LIBRARIES} Threads::Threads WavHelper AudioHelper kissfft FFTWrapper Util ParticleFilter nlohmann_json::nlohmann_json MapRenderer AudioCodec Config RingBuffer) 
# target_link_libraries(demo PRIVATE ${PORTAUDIO_LIBRARIES})
target_link_libraries(demo PRIVATE portaudio)
#target_link_libraries(demo PRIVATE ${SDL2_LIBRARIES} ${SDL2_TTF_LIBRARIES} )
target_link_libraries(demo PRIVATE SDL2 SDL2_ttf)
target_link_libraries(demo PRIVATE spdlog::spdlog)



# Copy gnuplot-iostream.h to the build directory
#configure_file(src/gnuplot-iostream.h ${CMAKE_CURRENT_BINARY_DIR}/gnuplot-iostream.h COPYONLY)
