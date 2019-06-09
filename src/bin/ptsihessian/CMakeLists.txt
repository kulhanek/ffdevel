# ==============================================================================
# FFDevel CMake File
# ==============================================================================

# program objects --------------------------------------------------------------
SET(PTSIHESS_SRC
        ptsihessian.f90
        )

# final build ------------------------------------------------------------------
ADD_EXECUTABLE(ptsihessian ${PTSIHESS_SRC})
ADD_DEPENDENCIES(ptsihessian ffdevel_shared)

TARGET_LINK_LIBRARIES(ptsihessian
        ffdevel
        ${PRMFILE_FLIB_NAME}
        ${SCIMAFIC_FLIB_NAME}
        ${LAPACK_LIB_NAME}
        ${BLAS_LIB_NAME}
        )

INSTALL(TARGETS
            ptsihessian
        DESTINATION
            bin
        )
