# define core ice model library
add_library(stim "")

if(GOTM_USE_STIM)
   project(gotm_stim VERSION 0.1.0 LANGUAGES Fortran)
   include(${CMAKE_CURRENT_LIST_DIR}/drivers/gotm/CMakeLists.txt)
   set(GET_GIT_INFO_SCRIPT "${gotm_SOURCE_DIR}/cmake/Modules/GetGitInfo.cmake")
   add_custom_target(stim_version
      ${CMAKE_COMMAND} -DINFILE=${CMAKE_CURRENT_LIST_DIR}/stim_version.F90.in -DOUTFILE=${CMAKE_CURRENT_BINARY_DIR}/src/stim_version.F90 -P ${GET_GIT_INFO_SCRIPT}
      WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}
      COMMENT "Retrieving description of last STIM commit..."
      SOURCES ${CMAKE_CURRENT_LIST_DIR}/stim_version.F90.in ${GET_GIT_INFO_SCRIPT}
      VERBATIM
   )
   add_dependencies(stim stim_version)
   set_property(TARGET stim_version PROPERTY FOLDER stim)
endif(GOTM_USE_STIM)

include(${CMAKE_CURRENT_LIST_DIR}/models/CMakeLists.txt)

set_property(TARGET stim PROPERTY FOLDER stim)

install(TARGETS stim
        EXPORT stimConfig
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
)
install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules/\${BUILD_TYPE}/
        DESTINATION include
)
