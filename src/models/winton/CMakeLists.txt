if(STIM_WINTON)
   target_sources(stim
       PRIVATE
       ${CMAKE_CURRENT_LIST_DIR}/ice_thm.F90
       ${CMAKE_CURRENT_LIST_DIR}/stim_winton.F90
   )
   target_compile_definitions(stim PUBLIC STIM_WINTON)
endif(STIM_WINTON)
