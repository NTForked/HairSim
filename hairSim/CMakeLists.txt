include( CMakeSourceFiles.txt )

# Eigen 3 library is required
find_package( Eigen3 REQUIRED )
include_directories( ${EIGEN3_INCLUDE_DIR} )

# MKL or LAPACK is required
if( NOT WETA STREQUAL "ON" )
  set( ENV{BLA_VENDOR} ${SUGGESTED_BLAS_LAPACK} )
  # This should also pick up BLAS
  find_package( LAPACK REQUIRED )
else( NOT WETA STREQUAL "ON" )
  find_package( MKL )
  if( MKL_FOUND)
      add_definitions(-DWITH_MKL)
      set ( LAPACK_INCLUDE_DIR ${MKL_INCLUDES} )
      set ( LAPACK_LIBRARIES ${MKL_LIBRARIES} )
  else( MKL_FOUND )
      set( ENV{BLA_VENDOR} ${SUGGESTED_BLAS_LAPACK} )
      find_package( LAPACK REQUIRED )
  endif( MKL_FOUND )
  include_directories( ${LAPACK_INCLUDE_DIR} )
endif( NOT WETA STREQUAL "ON" )

# Boost is required
# Note that prior to 1.43 there is a bug in boost caught by gcc-4.5 and up in c++0x
find_package( Boost 1.43.0 COMPONENTS serialization thread system REQUIRED )
include_directories( ${Boost_INCLUDE_DIRS} )

# OpenGL is required
find_package( OpenGL REQUIRED )
#include_directories( ${OPENGL_INCLUDE_DIR} )

# GLUT is required
find_package (GLUT REQUIRED glut)
if (GLUT_FOUND)
  include_directories (${GLUT_INCLUDE_DIR})
else (GLUT_FOUND)
  message (SEND_ERROR "Unable to locate GLUT")
endif (GLUT_FOUND)

# OpenMP is optional
find_package( OpenMP )
if( OPENMP_FOUND )
  set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}" )
endif( OPENMP_FOUND )

find_package (TCLAP REQUIRED)
include_directories (${TCLAP_INCLUDE_DIR})

# hairSim
add_executable( hairSimApp hairSimApp.cpp ${Headers} ${Sources} ${FORTRAN_SOURCES} )

# Mac OSX
target_link_libraries( hairSimApp bogus ${LAPACK_LIBRARIES} ${Boost_SYSTEM_LIBRARY} ${Boost_THREAD_LIBRARY} ${Boost_SERIALIZATION_LIBRARY} ${OPENGL_LIBRARIES} ${GLUT_glut_LIBRARY} )

# if Linux:
#target_link_libraries( hairSimApp bogus ${LAPACK_LIBRARIES} ${Boost_SYSTEM_LIBRARY} ${Boost_THREAD_LIBRARY} ${Boost_SERIALIZATION_LIBRARY} ${OPENGL_LIBRARIES} ${GLUT_glut_LIBRARY} )


