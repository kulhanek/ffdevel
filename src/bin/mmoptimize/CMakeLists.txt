# ==============================================================================
# FFDevel CMake File
# ==============================================================================

# program objects --------------------------------------------------------------
SET(MMOPTIMIZE_SRC
        mmoptimize.f90
        )

# final build ------------------------------------------------------------------
ADD_EXECUTABLE(mmoptimize ${MMOPTIMIZE_SRC})

TARGET_LINK_LIBRARIES(mmoptimize
        ffdevel
        ${PRMFILE_FLIB_NAME}
        ${SCIMAFIC_FLIB_NAME}
        ${LAPACK_LIB_NAME}
        ${BLAS_LIB_NAME}
        )

INSTALL(TARGETS
            mmoptimize
        DESTINATION
            bin
        )
