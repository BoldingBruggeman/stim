cmake_minimum_required(VERSION 3.0)

message("Building test_stim executable")
add_executable(test_stim EXCLUDE_FROM_ALL
               test_stim.F90
              )
target_link_libraries(test_stim PRIVATE stim)

set_property(TARGET test_stim PROPERTY FOLDER stim)


