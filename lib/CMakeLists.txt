# Define the library
add_library(WavHelper SHARED WavHelper/wavHelper.cpp ${CMAKE_CURRENT_SOURCE_DIR}/../src/main.h)
add_library(AudioHelper SHARED AudioHelper/audioHelper.cpp ${CMAKE_CURRENT_SOURCE_DIR}/../src/main.h)
add_library(kissfft SHARED kissfft/kiss_fft.c kissfft/kiss_fftr.c)
add_library(FFTWrapper SHARED FFTWrapper/fftWrapper.cpp ${CMAKE_CURRENT_SOURCE_DIR}/../src/main.h)
add_library(Util SHARED Util/util.cpp ${CMAKE_CURRENT_SOURCE_DIR}/../src/main.h)
add_library(ParticleFilter SHARED ParticleFilter/particleFilter.cpp  ${CMAKE_CURRENT_SOURCE_DIR}/../src/main.h ParticleFilter/particle.cpp ParticleFilter/Map/mapData.h ParticleFilter/Map/cell.cpp ParticleFilter/Map/wall.cpp ParticleFilter/Map/door.cpp ) 
add_library(MapRenderer SHARED MapRenderer/mapRenderer.cpp)

# Set the library output directory
set_target_properties(WavHelper PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)
set_target_properties(AudioHelper PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)
set_target_properties(FFTWrapper PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)
set_target_properties(Util PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)
set_target_properties(ParticleFilter PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)
set_target_properties(MapRenderer PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)
 
# Include the header files in the library target
target_include_directories(WavHelper PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/WavHelper ${CMAKE_CURRENT_SOURCE_DIR}/../src) 
target_include_directories(AudioHelper PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/AudioHelper ${CMAKE_CURRENT_SOURCE_DIR}/../src) 
target_include_directories(kissfft PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/kissfft)
target_include_directories(FFTWrapper PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/FFTWrapper ${CMAKE_CURRENT_SOURCE_DIR}/../src) 
target_include_directories(Util PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/Util ${CMAKE_CURRENT_SOURCE_DIR}/../src) 
target_include_directories(ParticleFilter PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/ParticleFilter ${CMAKE_CURRENT_SOURCE_DIR}/../src) 
target_include_directories(MapRenderer PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/MapRenderer ${CMAKE_CURRENT_SOURCE_DIR}/../src) 

#Linking libraries:
target_link_libraries(FFTWrapper PUBLIC kissfft)
target_link_libraries(AudioHelper PUBLIC Util)
target_link_libraries(WavHelper PUBLIC Util)
target_link_libraries(ParticleFilter PUBLIC Util nlohmann_json::nlohmann_json)
target_link_libraries(MapRenderer PUBLIC Util ParticleFilter ${SDL2_LIBRARIES} ${SDL2_TTF_LIBRARIES})