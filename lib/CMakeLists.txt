# Define the library
add_library(WavReader SHARED WavReader/wavReader.cpp)
add_library(AudioHelper SHARED AudioHelper/audioHelper.cpp ${CMAKE_CURRENT_SOURCE_DIR}/../src/main.h)
add_library(kissfft SHARED kissfft/kiss_fft.c kissfft/kiss_fftr.c)
add_library(FFTWrapper SHARED FFTWrapper/fftWrapper.cpp ${CMAKE_CURRENT_SOURCE_DIR}/../src/main.h)

# Set the library output directory
set_target_properties(WavReader PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)
set_target_properties(AudioHelper PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)
set_target_properties(FFTWrapper PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)
 
# Include the header files in the library target
target_include_directories(WavReader PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/WavReader) 
target_include_directories(AudioHelper PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/AudioHelper ${CMAKE_CURRENT_SOURCE_DIR}/../src) 
target_include_directories(kissfft PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/kissfft)
target_include_directories(FFTWrapper PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/FFTWrapper ${CMAKE_CURRENT_SOURCE_DIR}/../src) 

#Linking libraries:
target_link_libraries(FFTWrapper PUBLIC kissfft)