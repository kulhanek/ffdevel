# ==============================================================================
# FFDevel CMake File
# ==============================================================================

# program objects --------------------------------------------------------------
SET(DENSOVERLAP_SRC
        densoverlap.f90
        )

# final build ------------------------------------------------------------------
ADD_EXECUTABLE(densoverlap ${DENSOVERLAP_SRC})
ADD_DEPENDENCIES(densoverlap ffdevel_shared)

TARGET_LINK_LIBRARIES(densoverlap
        ffdevel
        ${PRMFILE_FLIB_NAME}
        ${SCIMAFIC_FLIB_NAME}
        ${LAPACK_LIB_NAME}
        ${BLAS_LIB_NAME}
        )

INSTALL(TARGETS
            densoverlap
        DESTINATION
            bin
        )
