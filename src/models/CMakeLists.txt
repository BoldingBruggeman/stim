# add definitions and sources to the core ice model library
target_sources(stim
    PRIVATE
    ${CMAKE_CURRENT_BINARY_DIR}/src/stim_version.F90
    ${CMAKE_CURRENT_LIST_DIR}/stim_variables.F90
    ${CMAKE_CURRENT_LIST_DIR}/stim_models.F90
)
set_property(TARGET stim PROPERTY Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules)
include(${CMAKE_CURRENT_LIST_DIR}/lebedev/CMakeLists.txt)
include(${CMAKE_CURRENT_LIST_DIR}/mylake/CMakeLists.txt)
include(${CMAKE_CURRENT_LIST_DIR}/winton/CMakeLists.txt)
#KBinclude(${CMAKE_CURRENT_LIST_DIR}/obsice/CMakeLists.txt)
include(${CMAKE_CURRENT_LIST_DIR}/basal_melt/CMakeLists.txt)
target_include_directories(stim PUBLIC
        $<INSTALL_INTERFACE:include>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/modules>
)

