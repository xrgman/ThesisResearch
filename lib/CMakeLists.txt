# Define the library
add_library(WavHelper SHARED WavHelper/wavHelper.cpp ${CMAKE_CURRENT_SOURCE_DIR}/../src/main.h)
add_library(AudioHelper SHARED AudioHelper/audioHelper.cpp ${CMAKE_CURRENT_SOURCE_DIR}/../src/main.h)
add_library(kissfft SHARED kissfft/kiss_fft.c kissfft/kiss_fftr.c)
add_library(FFTWrapper SHARED FFTWrapper/fftWrapper.cpp ${CMAKE_CURRENT_SOURCE_DIR}/../src/main.h)
add_library(Util SHARED Util/util.cpp ${CMAKE_CURRENT_SOURCE_DIR}/../src/main.h)

# Set the library output directory
set_target_properties(WavHelper PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)
set_target_properties(AudioHelper PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)
set_target_properties(FFTWrapper PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)
set_target_properties(Util PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)
 
# Include the header files in the library target
target_include_directories(WavHelper PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/WavHelper ${CMAKE_CURRENT_SOURCE_DIR}/../src) 
target_include_directories(AudioHelper PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/AudioHelper ${CMAKE_CURRENT_SOURCE_DIR}/../src) 
target_include_directories(kissfft PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/kissfft)
target_include_directories(FFTWrapper PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/FFTWrapper ${CMAKE_CURRENT_SOURCE_DIR}/../src) 
target_include_directories(Util PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/Util ${CMAKE_CURRENT_SOURCE_DIR}/../src) 

#Linking libraries:
target_link_libraries(FFTWrapper PUBLIC kissfft)
target_link_libraries(AudioHelper PUBLIC Util)