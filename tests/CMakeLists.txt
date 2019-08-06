set(ICE_BUILD_LIBRARIES_ONLY off)

if(NOT ICE_BUILD_LIBRARIES_ONLY)
message("Building test_stim executable")
add_executable(test_stim
               test_stim.F90
              )
target_link_libraries(test_stim PRIVATE stim)

install(TARGETS test_stim DESTINATION bin)
endif(NOT ICE_BUILD_LIBRARIES_ONLY)

