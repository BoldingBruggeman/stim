# define core ice model library
add_library(stim "")

if(GOTM_USE_STIM)
   project(gotm_stim VERSION 0.1.0 LANGUAGES Fortran)
   include(${CMAKE_CURRENT_LIST_DIR}/drivers/gotm/CMakeLists.txt)
endif(GOTM_USE_STIM)

include(${CMAKE_CURRENT_LIST_DIR}/models/CMakeLists.txt)

install(TARGETS stim
        EXPORT stimConfig
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
)
install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules/\${BUILD_TYPE}/
        DESTINATION include
)
