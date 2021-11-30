if(STIM_BASAL_MELT)
   target_sources(stim 
       PRIVATE
       ${CMAKE_CURRENT_LIST_DIR}/stim_basal_melt.F90
   )
   target_compile_definitions(stim PUBLIC STIM_BASAL_MELT)
endif(STIM_BASAL_MELT)
