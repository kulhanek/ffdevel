# ==============================================================================
# FFDevel CMake File
# ==============================================================================

# program objects --------------------------------------------------------------
SET(RSTGRADIENT_SRC
        rstgradient.f90
        )

# final build ------------------------------------------------------------------
ADD_EXECUTABLE(rstgradient ${RSTGRADIENT_SRC})
ADD_DEPENDENCIES(rstgradient ffdevel_shared)

TARGET_LINK_LIBRARIES(rstgradient
        ffdevel
        ${PRMFILE_FLIB_NAME}
        ${SCIMAFIC_FLIB_NAME}
        ${LAPACK_LIB_NAME}
        ${BLAS_LIB_NAME}
        )

INSTALL(TARGETS
            rstgradient
        DESTINATION
            bin
        )
