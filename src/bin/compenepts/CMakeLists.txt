# ==============================================================================
# FFDevel CMake File
# ==============================================================================

# program objects --------------------------------------------------------------
SET(COMP_ENE_SRC
        compenepts.f90
        )

# final build ------------------------------------------------------------------
ADD_EXECUTABLE(compenepts ${COMP_ENE_SRC})
ADD_DEPENDENCIES(compenepts ffdevel_shared)

TARGET_LINK_LIBRARIES(compenepts
        ffdevel
        ${PRMFILE_FLIB_NAME}
        ${SCIMAFIC_FLIB_NAME}
        ${LAPACK_LIB_NAME}
        ${BLAS_LIB_NAME}
        )

INSTALL(TARGETS
            compenepts
        DESTINATION
            bin
        )
