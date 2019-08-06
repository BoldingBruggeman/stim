set(ICE_BUILD_LIBRARIES_ONLY off)

if(NOT ICE_BUILD_LIBRARIES_ONLY)
message("Building test_ice executable")
add_executable(test_ice
               test_ice.F90
              )
target_link_libraries(test_ice PRIVATE ice)

install(TARGETS test_ice DESTINATION bin)
endif(NOT ICE_BUILD_LIBRARIES_ONLY)

