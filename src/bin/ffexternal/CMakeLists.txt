# ==============================================================================
# FFDevel CMake File
# ==============================================================================

# program objects --------------------------------------------------------------
SET(FFEXTERNAL_SRC
        ffexternal.f90
        )

# final build ------------------------------------------------------------------
ADD_EXECUTABLE(ffexternal ${FFEXTERNAL_SRC})

TARGET_LINK_LIBRARIES(ffexternal
        ffdevel
        ${PRMFILE_FLIB_NAME}
        ${SCIMAFIC_FLIB_NAME}
        ${LAPACK_LIB_NAME}
        ${BLAS_LIB_NAME}
        )

INSTALL(TARGETS
            ffexternal
        DESTINATION
            bin
        )
