cmake_minimum_required(VERSION 3.0)

if(GOTM_USE_ICE)
   project(gotm_stim VERSION 0.1.0 LANGUAGES Fortran)
   set(driver_taget "gotm_stim")
endif(GOTM_USE_ICE)

if(GETM_USE_ICE)
   project(getm_stim VERSION 0.1.0 LANGUAGES Fortran)
   set(driver_taget "getm_stim")
endif(GETM_USE_ICE)

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules)

# Selecting driver model

if(GOTM_USE_ICE)
add_library(gotm_stim
            ${CMAKE_CURRENT_BINARY_DIR}/stim_version.F90
            drivers/gotm_stim.F90
           )
set_property(TARGET gotm_stim PROPERTY Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules)
target_include_directories(gotm_stim PUBLIC
        $<INSTALL_INTERFACE:include>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/modules>
)
target_link_libraries(gotm_stim PRIVATE stim)
endif(GOTM_USE_ICE)

if(GETM_USE_ICE)
add_library(getm_stim
            drivers/getm_stim.F90
           )
set_property(TARGET getm_stim PROPERTY Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules)
target_include_directories(getm_stim PUBLIC
        $<INSTALL_INTERFACE:include>
        $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}>
        $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/modules>
)
target_link_libraries(getm_stim PRIVATE stim)
endif(GETM_USE_ICE)

if(STIM_LEBEDEV)
   set(SRC_LEBEDEV models/stim_lebedev.F90)
endif(STIM_LEBEDEV)
if(STIM_MYLAKE)
   set(SRC_MYLAKE models/stim_mylake.F90)
endif(STIM_MYLAKE)
if(STIM_WINTON)
   set(SRC_WINTON models/ice_thm.F90 models/stim_winton.F90)
endif(STIM_WINTON)

add_library(stim
            models/stim_variables.F90
            models/stim_models.F90
            ${SRC_LEBEDEV}
            ${SRC_MYLAKE}
            ${SRC_WINTON}
           )
if(STIM_LEBEDEV)
   target_compile_definitions(stim PUBLIC STIM_LEBEDEV)
endif(STIM_LEBEDEV)
if(STIM_MYLAKE)
   target_compile_definitions(stim PUBLIC STIM_MYLAKE)
endif(STIM_MYLAKE)
if(STIM_WINTON)
   target_compile_definitions(stim PUBLIC STIM_WINTON)
endif(STIM_WINTON)
set_property(TARGET stim PROPERTY Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules)
target_include_directories(stim PUBLIC
        $<INSTALL_INTERFACE:include>
        $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/modules>
)

install(TARGETS gotm_stim stim
        EXPORT stimConfig
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
)
install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules/\${BUILD_TYPE}/
        DESTINATION include
)
       
