cmake_minimum_required(VERSION 3.0)

# Create main project.
project(stim VERSION 0.1.0 LANGUAGES Fortran)

# Use solution folders in IDEs
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

# Create placeholder file with version information
# used when STIM_EMBED_VERSION is off.
set(GIT_COMMIT_ID ${VERSION})
set(GIT_BRANCH_NAME unknown)
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/src/stim_version.F90.in" ${CMAKE_CURRENT_BINARY_DIR}/src/stim_version.F90)

# Use use position-independent code (-fPIC) everywhere if building shared libraries
if(BUILD_SHARED_LIBS)
  set(CMAKE_POSITION_INDEPENDENT_CODE ON)
endif()

#[[---
option(STIM_EMBED_VERSION "Include Git version info" "ON")
if(STIM_EMBED_VERSION)
   add_custom_target(version ALL
      ${CMAKE_COMMAND} -DINFILE=${PROJECT_SOURCE_DIR}/src/stim_version.F90.in -DOUTFILE=${CMAKE_CURRENT_BINARY_DIR}/src/stim_version.F90 -P "${PROJECT_SOURCE_DIR}/cmake/Modules/GetGitInfo.cmake"
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
      COMMENT "Retrieving description of last ICE commit..."
      SOURCES ${PROJECT_SOURCE_DIR}/src/stim_version.F90.in "${PROJECT_SOURCE_DIR}/cmake/Modules/GetGitInfo.cmake"
      VERBATIM
  )
endif(STIM_EMBED_VERSION)
---]]

option(STIM_LEBEDEV "Lebedev ICE model" ON)
option(STIM_MYLAKE "MyLake ICE model" ON)
option(STIM_WINTON "Winton ICE model" ON)

add_subdirectory(src)
add_subdirectory(tests)

install(EXPORT stimConfig DESTINATION cmake)
