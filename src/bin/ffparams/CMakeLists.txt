# ==============================================================================
# FFDevel CMake File
# ==============================================================================

# program objects --------------------------------------------------------------
SET(FFPARAMS_SRC
        ffparams.f90
        )

# final build ------------------------------------------------------------------
ADD_EXECUTABLE(ffparams ${FFPARAMS_SRC})
ADD_DEPENDENCIES(ffparams ffdevel_shared)

TARGET_LINK_LIBRARIES(ffparams
        ffdevel
        ${PRMFILE_FLIB_NAME}
        ${SCIMAFIC_FLIB_NAME}
        ${LAPACK_LIB_NAME}
        ${BLAS_LIB_NAME}
        )

INSTALL(TARGETS
            ffparams
        DESTINATION
            bin
        )
