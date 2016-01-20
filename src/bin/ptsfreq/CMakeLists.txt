# ==============================================================================
# FFDevel CMake File
# ==============================================================================

# program objects --------------------------------------------------------------
SET(PTSFREQ_SRC
        ptsfreq.f90
        )

# final build ------------------------------------------------------------------
ADD_EXECUTABLE(ptsfreq ${PTSFREQ_SRC})

TARGET_LINK_LIBRARIES(ptsfreq
        ffdevel
        ${PRMFILE_FLIB_NAME}
        ${SCIMAFIC_FLIB_NAME}
        ${LAPACK_LIB_NAME}
        ${BLAS_LIB_NAME}
        )

INSTALL(TARGETS
            ptsfreq
        DESTINATION
            bin
        )
