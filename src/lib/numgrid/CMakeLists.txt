# ==============================================================================
# NumGrid CMake File
# ==============================================================================

SET(NUMGRID_VERSION 1.1.0)

INCLUDE_DIRECTORIES(lebedev SYSTEM)

# objects in library -----------------------------------------------------------
SET(NUMGRID_SRC
        AtomGrid.cpp
        becke_partitioning.cpp
        bragg.cpp
        grid_radial.cpp
        version.cpp
        numgrid.f90
        lebedev/sphere_lebedev_rule.cpp
        )

# create static library --------------------------------------------------------
IF(LIBS_STATIC)
    ADD_LIBRARY(numgrid_static STATIC ${NUMGRID_SRC})
    SET_TARGET_PROPERTIES(numgrid_static PROPERTIES
                            OUTPUT_NAME numgrid
                            CLEAN_DIRECT_OUTPUT 1)
ENDIF(LIBS_STATIC)

# create shared library --------------------------------------------------------
IF(LIBS_SHARED)
    ADD_LIBRARY(numgrid_shared SHARED ${NUMGRID_SRC})
    SET_TARGET_PROPERTIES(numgrid_shared PROPERTIES
                            OUTPUT_NAME numgrid
                            CLEAN_DIRECT_OUTPUT 1
                            VERSION ${NUMGRID_VERSION})
    INSTALL(TARGETS
                numgrid_shared
            LIBRARY DESTINATION
                lib
            )
ENDIF(LIBS_SHARED)
