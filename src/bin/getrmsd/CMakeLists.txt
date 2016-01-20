# ==============================================================================
# FFDevel CMake File
# ==============================================================================

# program objects --------------------------------------------------------------
SET(GETRMSD_SRC
        getrmsd.f90
        )

# final build ------------------------------------------------------------------
ADD_EXECUTABLE(getrmsd ${GETRMSD_SRC})
ADD_DEPENDENCIES(getrmsd ffdevel_shared)

TARGET_LINK_LIBRARIES(getrmsd
        ffdevel
        ${PRMFILE_FLIB_NAME}
        ${SCIMAFIC_FLIB_NAME}
        ${LAPACK_LIB_NAME}
        ${BLAS_LIB_NAME}
        )

INSTALL(TARGETS
            getrmsd
        DESTINATION
            bin
        )
