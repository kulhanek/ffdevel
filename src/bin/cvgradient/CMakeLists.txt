# ==============================================================================
# FFDevel CMake File
# ==============================================================================

# program objects --------------------------------------------------------------
SET(CVGRADIENT_SRC
        cvgradient.f90
        )

# final build ------------------------------------------------------------------
ADD_EXECUTABLE(cvgradient ${CVGRADIENT_SRC})
ADD_DEPENDENCIES(cvgradient ffdevel_shared)

TARGET_LINK_LIBRARIES(cvgradient
        ffdevel
        ${PRMFILE_FLIB_NAME}
        ${SCIMAFIC_FLIB_NAME}
        ${LAPACK_LIB_NAME}
        ${BLAS_LIB_NAME}
        )

INSTALL(TARGETS
            cvgradient
        DESTINATION
            bin
        )
