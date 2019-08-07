target_sources(stim
    PRIVATE
    ${CMAKE_CURRENT_BINARY_DIR}/src/stim_version.F90
    ${CMAKE_CURRENT_LIST_DIR}/stim_variables.F90
    ${CMAKE_CURRENT_LIST_DIR}/stim_models.F90
)
include(${CMAKE_CURRENT_LIST_DIR}/lebedev/CMakeLists.txt)
include(${CMAKE_CURRENT_LIST_DIR}/mylake/CMakeLists.txt)
include(${CMAKE_CURRENT_LIST_DIR}/winton/CMakeLists.txt)

