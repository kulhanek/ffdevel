# ==============================================================================
# FFDevel CMake File
# ==============================================================================

# program objects --------------------------------------------------------------
SET(FFHESSIAN_SRC
        ffhessian.f90
        )

# final build ------------------------------------------------------------------
ADD_EXECUTABLE(ffhessian ${FFHESSIAN_SRC})

TARGET_LINK_LIBRARIES(ffhessian
        ffdevel
        ${PRMFILE_FLIB_NAME}
        ${SCIMAFIC_FLIB_NAME}
        ${LAPACK_LIB_NAME}
        ${BLAS_LIB_NAME}
        )

INSTALL(TARGETS
            ffhessian
        DESTINATION
            bin
        )
