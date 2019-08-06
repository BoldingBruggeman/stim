cmake_minimum_required(VERSION 3.10)

# Create main project.
project(ice VERSION 0.1.0 LANGUAGES Fortran)

# Create placeholder file with version information
# used when ICE_EMBED_VERSION is off.
set(GIT_COMMIT_ID ${VERSION})
set(GIT_BRANCH_NAME unknown)
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/src/ice_version.F90.in" ${CMAKE_CURRENT_BINARY_DIR}/src/ice_version.F90)

#[[---
if(ICE_EMBED_VERSION)
   add_custom_target(version ALL
      ${CMAKE_COMMAND} -DINFILE=${PROJECT_SOURCE_DIR}/src/ice_version.F90.in -DOUTFILE=${CMAKE_CURRENT_BINARY_DIR}/src/ice_version.F90 -P "${PROJECT_SOURCE_DIR}/cmake/Modules/GetGitInfo.cmake"
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
      COMMENT "Retrieving description of last ICE commit..."
      SOURCES ${PROJECT_SOURCE_DIR}/src/ice_version.F90.in "${PROJECT_SOURCE_DIR}/cmake/Modules/GetGitInfo.cmake"
      VERBATIM
  )
endif(ICE_EMBED_VERSION)
---]]

option(ICE_LEBEDEV "Lebedev ICE model" ON)
option(ICE_MYLAKE "MyLake ICE model" ON)
option(ICE_WINTON "Winton ICE model" ON)

add_subdirectory(src)
add_subdirectory(tests EXCLUDE_FROM_DEFAULT)
