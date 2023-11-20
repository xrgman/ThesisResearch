# Define the library
add_library(WavReader SHARED WavReader/wavReader.cpp)
add_library(AudioHelper SHARED AudioHelper/audioHelper.cpp)

# Set the library output directory
set_target_properties(WavReader PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)
set_target_properties(AudioHelper PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)
 
# Include the header files in the library target
target_include_directories(WavReader PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/WavReader) 
target_include_directories(AudioHelper PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/AudioHelper) 