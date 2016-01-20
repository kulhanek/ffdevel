# ==============================================================================
# FFDevel CMake File
# ==============================================================================

# program objects --------------------------------------------------------------
SET(RMSDFIT_SRC
        rmsdfit.f90
        )

# final build ------------------------------------------------------------------
ADD_EXECUTABLE(rmsdfit ${RMSDFIT_SRC})
ADD_DEPENDENCIES(rmsdfit ffdevel_shared)

TARGET_LINK_LIBRARIES(rmsdfit
        ffdevel
        ${PRMFILE_FLIB_NAME}
        ${SCIMAFIC_FLIB_NAME}
        ${LAPACK_LIB_NAME}
        ${BLAS_LIB_NAME}
        )

INSTALL(TARGETS
            rmsdfit
        DESTINATION
            bin
        )
