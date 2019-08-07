add_library(gotm_stim "")
target_sources(gotm_stim
    PRIVATE
    ${CMAKE_CURRENT_LIST_DIR}/gotm_stim.F90
)
target_include_directories(gotm_stim
    PUBLIC
        $<INSTALL_INTERFACE:include>
        $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/modules>
)
target_link_libraries(gotm_stim PRIVATE stim)

