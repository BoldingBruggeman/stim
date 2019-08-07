# define core ice model library
add_library(stim "")

if(GOTM_USE_STIM)
   project(gotm_stim VERSION 0.1.0 LANGUAGES Fortran)
   include(${CMAKE_CURRENT_LIST_DIR}/drivers/gotm/CMakeLists.txt)
   set_property(TARGET gotm_stim stim PROPERTY Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules)
   target_link_libraries(gotm_stim PRIVATE stim)
endif(GOTM_USE_STIM)

# add definitions and sources to the core ice model library
target_include_directories(stim PUBLIC
        $<INSTALL_INTERFACE:include>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/modules>
)
include(${CMAKE_CURRENT_LIST_DIR}/models/CMakeLists.txt)

install(TARGETS stim
        EXPORT stimConfig
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
)
install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules/\${BUILD_TYPE}/
        DESTINATION include
)
