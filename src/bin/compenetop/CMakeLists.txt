# ==============================================================================
# FFDevel CMake File
# ==============================================================================

# program objects --------------------------------------------------------------
SET(COMP_ENE_SRC
        compenetop.f90
        )

# final build ------------------------------------------------------------------
ADD_EXECUTABLE(compenetop ${COMP_ENE_SRC})

TARGET_LINK_LIBRARIES(compenetop
        ffdevel
        ${PRMFILE_FLIB_NAME}
        ${SCIMAFIC_FLIB_NAME}
        ${LAPACK_LIB_NAME}
        ${BLAS_LIB_NAME}
        )

INSTALL(TARGETS
            compenetop
        DESTINATION
            bin
        )
