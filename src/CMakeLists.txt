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

if(ICE_LEBEDEV)
   set(SRC_LEBEDEV models/ice_lebedev.F90)
endif(ICE_LEBEDEV)
if(ICE_MYLAKE)
   set(SRC_MYLAKE models/ice_mylake.F90)
endif(ICE_MYLAKE)
if(ICE_WINTON)
   set(SRC_WINTON models/ice_thm.F90 models/ice_winton.F90)
endif(ICE_WINTON)

add_library(ice)
target_sources(ice
            PRIVATE
            models/ice_variables.F90
            models/ice_models.F90
            ${SRC_LEBEDEV}
            ${SRC_MYLAKE}
            ${SRC_WINTON}
           )
if(ICE_LEBEDEV)
   target_compile_definitions(ice PUBLIC ICE_LEBEDEV)
endif(ICE_LEBEDEV)
if(ICE_MYLAKE)
   target_compile_definitions(ice PUBLIC ICE_MYLAKE)
endif(ICE_MYLAKE)
if(ICE_WINTON)
   target_compile_definitions(ice PUBLIC ICE_WINTON)
endif(ICE_WINTON)
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
