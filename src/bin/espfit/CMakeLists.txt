# ==============================================================================
# FFDevel CMake File
# ==============================================================================

# program objects --------------------------------------------------------------
SET(FFOPTIMIZE_SRC
        ffoptimize.f90
        )

# final build ------------------------------------------------------------------
ADD_EXECUTABLE(ffoptimize ${FFOPTIMIZE_SRC})

TARGET_LINK_LIBRARIES(ffoptimize
        ffdevel
        ${PRMFILE_FLIB_NAME}
        ${SCIMAFIC_FLIB_NAME}
        ${LAPACK_LIB_NAME}
        ${BLAS_LIB_NAME}
        )

INSTALL(TARGETS
            ffoptimize
        DESTINATION
            bin
        )
