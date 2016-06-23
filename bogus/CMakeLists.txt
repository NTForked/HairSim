
set(CMAKE_MACOSX_RPATH 1)

execute_process( COMMAND ../updateCMakeSources.sh WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} )
file(STRINGS CMakeSources.txt SRC_LIST)

find_package( Eigen3 REQUIRED )
include_directories( ${EIGEN3_INCLUDE_DIR} )

add_definitions( -DBOGUS_INSTANTIATE_3D_SOC)

if( WITH_DYNAMIC )
        add_definitions( -DBOGUS_INSTANTIATE_DYNAMIC_SOC )
endif()

#SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -S --verbose-asm " )

add_library( bogus SHARED ${SRC_LIST})
target_link_libraries(bogus ${EXTERNAL_LIBRARIES} )

FOREACH(item ${SRC_LIST})
  STRING(REGEX MATCH "^Interfaces/.+\\.hpp$" item ${item})
  IF(item)
	LIST(APPEND INTERFACES ${item})
  ENDIF(item)
ENDFOREACH(item ${SRC_LIST})

INSTALL( DIRECTORY Core DESTINATION include/bogus )
INSTALL( DIRECTORY Extra DESTINATION include/bogus )
INSTALL( FILES ${INTERFACES} DESTINATION include/bogus/Interfaces )
INSTALL( TARGETS bogus
  LIBRARY DESTINATION lib
  RUNTIME DESTINATION lib
  ARCHIVE DESTINATION lib )

