if(STIM_LEBEDEV)
   target_sources(stim 
       PRIVATE
       ${CMAKE_CURRENT_LIST_DIR}/stim_lebedev.F90
   )
   target_compile_definitions(stim PUBLIC STIM_LEBEDEV)
endif(STIM_LEBEDEV)
