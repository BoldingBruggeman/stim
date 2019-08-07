add_executable(test_stim EXCLUDE_FROM_ALL
               ${CMAKE_CURRENT_LIST_DIR}/test_stim.F90
              )
target_link_libraries(test_stim PRIVATE stim)
set_property(TARGET test_stim PROPERTY FOLDER stim)
