# ==============================================================================
# FFDevel CMake File
# ==============================================================================

# program objects --------------------------------------------------------------
SET(GETD3_SRC
        getd3.f90
        )

# final build ------------------------------------------------------------------
ADD_EXECUTABLE(getd3 ${GETD3_SRC})
ADD_DEPENDENCIES(getd3 ffdevel_shared)

TARGET_LINK_LIBRARIES(getd3
        ffdevel
        ${PRMFILE_FLIB_NAME}
        ${SCIMAFIC_FLIB_NAME}
        ${LAPACK_LIB_NAME}
        ${BLAS_LIB_NAME}
        )

INSTALL(TARGETS
            getd3
        DESTINATION
            bin
        )
