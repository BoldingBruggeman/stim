#add_library(ice OBJECT
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules)

# Selecting driver model

if(GOTM_USE_ICE)
add_library(gotm_ice)
target_sources(gotm_ice
            PRIVATE
            ${CMAKE_CURRENT_BINARY_DIR}/ice_version.F90
            drivers/gotm_ice.F90
           )
target_include_directories(gotm_ice PUBLIC
        $<INSTALL_INTERFACE:include>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/modules>
)
target_link_libraries(gotm_ice PRIVATE ice)
endif(GOTM_USE_ICE)

if(GETM_USE_ICE)
add_library(getm_ice)
target_sources(getm_ice
            PRIVATE
            drivers/getm_ice.F90
           )
target_include_directories(getm_ice PUBLIC
        $<INSTALL_INTERFACE:include>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/modules>
)
target_link_libraries(getm_ice PRIVATE ice)
endif(GETM_USE_ICE)

if(STIM_LEBEDEV)
   set(SRC_LEBEDEV models/ice_lebedev.F90)
endif(STIM_LEBEDEV)
if(STIM_MYLAKE)
   set(SRC_MYLAKE models/ice_mylake.F90)
endif(STIM_MYLAKE)
if(STIM_WINTON)
   set(SRC_WINTON models/ice_thm.F90 models/ice_winton.F90)
endif(STIM_WINTON)

add_library(ice)
target_sources(ice
            PRIVATE
            models/ice_variables.F90
            models/ice_models.F90
            ${SRC_LEBEDEV}
            ${SRC_MYLAKE}
            ${SRC_WINTON}
           )
if(STIM_LEBEDEV)
   target_compile_definitions(ice PUBLIC STIM_LEBEDEV)
endif(STIM_LEBEDEV)
if(STIM_MYLAKE)
   target_compile_definitions(ice PUBLIC STIM_MYLAKE)
endif(STIM_MYLAKE)
if(STIM_WINTON)
   target_compile_definitions(ice PUBLIC STIM_WINTON)
endif(STIM_WINTON)
target_include_directories(ice PUBLIC
        $<INSTALL_INTERFACE:include>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/modules>
)

install(FILES
#KB        ${CMAKE_CURRENT_BINARY_DIR}/modules/ice.mod
#KB        ${CMAKE_CURRENT_BINARY_DIR}/modules/ice_models.mod
#KB        ${CMAKE_CURRENT_BINARY_DIR}/modules/ice_variables.mod
        DESTINATION include)
install(TARGETS ice DESTINATION lib)
