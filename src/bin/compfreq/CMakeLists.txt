# ==============================================================================
# FFDevel CMake File
# ==============================================================================

# program objects --------------------------------------------------------------
SET(COMP_FREQ_SRC
        compfreq.f90
        )

# final build ------------------------------------------------------------------
ADD_EXECUTABLE(compfreq ${COMP_FREQ_SRC})

TARGET_LINK_LIBRARIES(compfreq
        ffdevel
        ${PRMFILE_FLIB_NAME}
        ${SCIMAFIC_FLIB_NAME}
        ${LAPACK_LIB_NAME}
        ${BLAS_LIB_NAME}
        )

INSTALL(TARGETS
            compfreq
        DESTINATION
            bin
        )
