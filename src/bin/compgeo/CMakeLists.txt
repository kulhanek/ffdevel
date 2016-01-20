# ==============================================================================
# FFDevel CMake File
# ==============================================================================

# program objects --------------------------------------------------------------
SET(COMP_GEO_SRC
        compgeo.f90
        )

# final build ------------------------------------------------------------------
ADD_EXECUTABLE(compgeo ${COMP_GEO_SRC})

TARGET_LINK_LIBRARIES(compgeo
        ffdevel
        ${PRMFILE_FLIB_NAME}
        ${SCIMAFIC_FLIB_NAME}
        ${LAPACK_LIB_NAME}
        ${BLAS_LIB_NAME}
        )

INSTALL(TARGETS
            compgeo
        DESTINATION
            bin
        )
