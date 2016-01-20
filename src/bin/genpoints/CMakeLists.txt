# ==============================================================================
# FFDevel CMake File
# ==============================================================================

# program objects --------------------------------------------------------------
SET(GENPOINTS_SRC
        genpoints.f90
        genpoints_dat.f90
        genpoints_control.f90
        genpoints_utils.f90
        genpoints_math.f90
        )

# final build ------------------------------------------------------------------
ADD_EXECUTABLE(genpoints ${GENPOINTS_SRC})

TARGET_LINK_LIBRARIES(genpoints
        ffdevel
        ${PRMFILE_FLIB_NAME}
        ${SCIMAFIC_FLIB_NAME}
        ${LAPACK_LIB_NAME}
        ${BLAS_LIB_NAME}
        )

INSTALL(TARGETS
            genpoints
        DESTINATION
            bin
        )
