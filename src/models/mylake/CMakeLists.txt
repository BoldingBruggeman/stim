if(STIM_MYLAKE)
   target_sources(stim 
       PRIVATE
       ${CMAKE_CURRENT_LIST_DIR}/stim_mylake.F90
   )
   target_compile_definitions(stim PUBLIC STIM_MYLAKE)
endif(STIM_MYLAKE)

