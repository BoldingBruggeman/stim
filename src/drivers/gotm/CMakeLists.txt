add_library(gotm_stim ${CMAKE_CURRENT_LIST_DIR}/gotm_stim.F90)
set_property(TARGET gotm_stim PROPERTY Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules)
set_property(TARGET gotm_stim PROPERTY FOLDER stim)
target_link_libraries(gotm_stim PRIVATE stim gotm::config)

