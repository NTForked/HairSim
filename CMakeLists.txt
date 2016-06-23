cmake_minimum_required( VERSION 2.8 )

project( hairSim C CXX Fortran )

###############################################################################
### Configure and check the build setting options

# If the user provided a CMAKE_BUILD_TYPE, create a cache variable
if( NOT CMAKE_BUILD_TYPE STREQUAL "" )
  set( CMAKE_BUILD_TYPE ${CMAKE_BUILD_TYPE} CACHE STRING "Choose the type of
       build, options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug
       Release RelWithDebInfo MinSizeRel." FORCE )

# Otherwise default to debug
else( NOT CMAKE_BUILD_TYPE STREQUAL "" )
  set( CMAKE_BUILD_TYPE Debug CACHE STRING "Choose the type of
       build, options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug
       Release RelWithDebInfo MinSizeRel." FORCE )
endif( NOT CMAKE_BUILD_TYPE STREQUAL "" )

# Enforce the use of a valid build type
if( NOT CMAKE_BUILD_TYPE STREQUAL "None" AND NOT CMAKE_BUILD_TYPE STREQUAL "Debug" )
  if( NOT CMAKE_BUILD_TYPE STREQUAL "Release" AND NOT CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo" )
    if( NOT CMAKE_BUILD_TYPE STREQUAL "MinSizeRel" )
      message (SEND_ERROR "Choose the type of build, options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug Release RelWithDebInfo MinSizeRel.")
    endif( NOT CMAKE_BUILD_TYPE STREQUAL "MinSizeRel" )
  endif( NOT CMAKE_BUILD_TYPE STREQUAL "Release" AND NOT CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo" )
endif( NOT CMAKE_BUILD_TYPE STREQUAL "None" AND NOT CMAKE_BUILD_TYPE STREQUAL "Debug" )

###############################################################################
### Allow the user to set a preferred BLAS and LAPACK implementation

set( SUGGESTED_BLAS_LAPACK "" CACHE STRING "The desired LAPACK implementation (e.g. Generic, Apple, ACML, Intel10_64lp, Intel10_64lp_seq, ..." )

###############################################################################
### Set flags used across all builds (TODO: Check compiler types here!)

# Append an extra space if flags are already set
if( CMAKE_CXX_FLAGS )
  set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} " )
endif( CMAKE_CXX_FLAGS )

# Enable aggressive warnings
set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wno-unused -lboost_serialization -lboost_system -Wno-unused-parameter -fopenmp" )

# Enable c++0x support
set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x" )

set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC" )

# give up after one error:
set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wfatal-errors" )

###############################################################################
### Set preprocessor defines used in all builds

add_definitions( -DEIGEN_DONT_PARALLELIZE )

###############################################################################
### Set preprocessor defines used in debug builds

if( CMAKE_BUILD_TYPE STREQUAL "Debug" )
  add_definitions( -DEIGEN_DONT_ALIGN ) # Also does -DEIGEN_DONT_VECTORIZE
endif( CMAKE_BUILD_TYPE STREQUAL "Debug" )

###############################################################################
### Set paths, etc

# Create a default location for custom cmake modules
set( CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake )

# Allow projects to reference source relative to the root of the project
include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )

###############################################################################
### Add libraries, plugins

add_subdirectory( bogus )
add_subdirectory( hairSim )
