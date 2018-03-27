CMAKE_MINIMUM_REQUIRED(VERSION 2.8.12)

set (CMAKE_CXX_STANDARD 11)

# Set the name of the project and target:
SET(TARGET "pitdealii")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic")


PROJECT(${TARGET})

set(PROJECT_PATH "${CMAKE_CURRENT_SOURCE_DIR}")
set(INCLUDE_PATH "${PROJECT_PATH}/src")
set(SRC_PATH "${PROJECT_PATH}/src")
set(TEST_SRC_PATH "${PROJECT_PATH}/test")

set(LIB_PATH  "${PROJECT_BINARY_DIR}/lib")
set(LIB_NAME  "${PROJECT_NAME}")
set(BIN_PATH  "${PROJECT_BINARY_DIR}/bin")
set(MAIN_NAME "${PROJECT_NAME}_main")
set(REFERENCE_NAME "${PROJECT_NAME}_reference")
set(TEST_PATH "${PROJECT_BINARY_DIR}/test")
set(TEST_NAME "test_${PROJECT_NAME}")

SET(MAIN_SRC ${SRC_PATH}/${TARGET}.cc)

SET(REFERENCE_SRC ${SRC_PATH}/reference.cc)

SET(SRC
  ${SRC_PATH}/BraidFuncs.cc
  ${SRC_PATH}/Utilities.cc
  )

SET(TEST_SRC
  ${TEST_SRC_PATH}/test_braid.cc
  )

OPTION(USE_MPI "Build MPI version" ON)

include_directories(AFTER ${INCLUDE_PATH})

if(USE_MPI)
find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})
add_definitions(-DUSE_MPI)
endif(USE_MPI)

FIND_PACKAGE(deal.II 9.0.0 QUIET
  HINTS ${deal.II_DIR} ${DEAL_II_DIR} $ENV{HOME}/Sources/builddealii ../../../builddealii ../ ../../ $ENV{DEAL_II_DIR}
  )
IF(NOT ${deal.II_FOUND})
  MESSAGE(FATAL_ERROR "\n"
    "*** Could not locate a (sufficiently recent) version of deal.II. ***\n\n"
    "You may want to either pass a flag -DDEAL_II_DIR=/path/to/deal.II to cmake\n"
    "or set an environment variable \"DEAL_II_DIR\" that contains this path."
    )
ENDIF()

DEAL_II_INITIALIZE_CACHED_VARIABLES()

# Find braid details
set(BRAID_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../braid")
# Include the braid paths and libraries
include_directories(${BRAID_DIR}/)


##############################
# Finally start building stuff

# First, build the "library" that consists of all the source files
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${LIB_PATH})
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${LIB_PATH})
add_library(${LIB_NAME} ${SRC})


# Next build the main function
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${BIN_PATH})
add_executable(${MAIN_NAME} ${MAIN_SRC})
target_link_libraries(${MAIN_NAME} ${LIB_NAME}) # link my library
target_link_libraries(${MAIN_NAME} ${BRAID_DIR}/libbraid.a) # link braid
DEAL_II_SETUP_TARGET(${MAIN_NAME}) # Link dealii
if(USE_MPI)
set_target_properties(${MAIN_NAME} PROPERTIES
  LINK_FLAGS "${MPI_LINK_FLAGS}") # Link MPI
set_target_properties(${MAIN_NAME} PROPERTIES
  COMPILE_FLAGS "${MPI_COMPILE_FLAGS}") # Use MPI compile flags
endif(USE_MPI)

# Build the reference
add_executable(${REFERENCE_NAME} ${REFERENCE_SRC})
target_link_libraries(${REFERENCE_NAME} ${LIB_NAME}) # link my library
DEAL_II_SETUP_TARGET(${REFERENCE_NAME}) # Link dealii
if(USE_MPI)
set_target_properties(${REFERENCE_NAME} PROPERTIES
  LINK_FLAGS "${MPI_LINK_FLAGS}")
set_target_properties(${REFERENCE_NAME} PROPERTIES
  COMPILE_FLAGS "${MPI_COMPILE_FLAGS}") # Use MPI compile flags
endif(USE_MPI)

# Lastly build the tests
enable_testing()
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${TEST_PATH})
add_subdirectory(test)

