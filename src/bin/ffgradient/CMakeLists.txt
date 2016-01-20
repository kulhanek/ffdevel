# ==============================================================================
# FFDevel CMake File
# ==============================================================================

# program objects --------------------------------------------------------------
SET(FFGRADIENT_SRC
        ffgradient.f90
        )

# final build ------------------------------------------------------------------
ADD_EXECUTABLE(ffgradient ${FFGRADIENT_SRC})
ADD_DEPENDENCIES(ffgradient ffdevel_shared)

TARGET_LINK_LIBRARIES(ffgradient
        ffdevel
        ${PRMFILE_FLIB_NAME}
        ${SCIMAFIC_FLIB_NAME}
        ${LAPACK_LIB_NAME}
        ${BLAS_LIB_NAME}
        )

INSTALL(TARGETS
            ffgradient
        DESTINATION
            bin
        )
