# ==============================================================================
# FFDevel CMake File
# ==============================================================================

# program objects --------------------------------------------------------------
SET(FFENERGY_SRC
        ffenergy.f90
        )

# final build ------------------------------------------------------------------
ADD_EXECUTABLE(ffenergy ${FFENERGY_SRC})
ADD_DEPENDENCIES(ffenergy ffdevel_shared)

TARGET_LINK_LIBRARIES(ffenergy
        ffdevel
        ${PRMFILE_FLIB_NAME}
        ${SCIMAFIC_FLIB_NAME}
        ${LAPACK_LIB_NAME}
        ${BLAS_LIB_NAME}
        )

INSTALL(TARGETS
            ffenergy
        DESTINATION
            bin
        )
