if(STIM_FLATO)
   target_sources(stim
       PRIVATE
       ${CMAKE_CURRENT_LIST_DIR}/ice_thm.F90
       ${CMAKE_CURRENT_LIST_DIR}/stim_flato.F90
   )
   target_compile_definitions(stim PUBLIC STIM_FLATO)
endif(STIM_FLATO)
