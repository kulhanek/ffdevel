# ==============================================================================
# FFDevel CMake File
# ==============================================================================

# program objects --------------------------------------------------------------
SET(VDWOPTIMIZE_SRC
        vdwoptimize.f90
        )

# final build ------------------------------------------------------------------
ADD_EXECUTABLE(vdwoptimize ${VDWOPTIMIZE_SRC})
ADD_DEPENDENCIES(vdwoptimize ffdevel_shared)

TARGET_LINK_LIBRARIES(vdwoptimize
        ffdevel
        ${PRMFILE_FLIB_NAME}
        ${SCIMAFIC_FLIB_NAME}
        ${LAPACK_LIB_NAME}
        ${BLAS_LIB_NAME}
        )

INSTALL(TARGETS
            vdwoptimize
        DESTINATION
            bin
        )
