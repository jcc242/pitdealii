# Set the name of the project and target:
SET(TARGET "pitdealii")

SET(TARGET_SRC
  ${TARGET}.cc
  )

# Usually, you will not need to modify anything beyond this point...

CMAKE_MINIMUM_REQUIRED(VERSION 2.8.12)

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

set(BRAID_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../braid")

DEAL_II_INITIALIZE_CACHED_VARIABLES()
PROJECT(${TARGET})
# DEAL_II_INVOKE_AUTOPILOT()
add_executable(${TARGET} ${TARGET_SRC})

DEAL_II_SETUP_TARGET(${TARGET})

# Include the braid paths and libraries
include_directories(${BRAID_DIR}/)
target_link_libraries(${TARGET} ${BRAID_DIR}/libbraid.a)
