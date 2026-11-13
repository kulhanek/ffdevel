# ==============================================================================
# FFDevel CMake File
# ==============================================================================

# program objects --------------------------------------------------------------
SET(RESP_SRC
    presp.f90
    )

# final build ------------------------------------------------------------------
ADD_EXECUTABLE(presp ${RESP_SRC})
TARGET_LINK_LIBRARIES(presp ${LAPACK_LIB_NAME} ${BLAS_LIB_NAME})

# install -----------------------------------------------------------------------
INSTALL(TARGETS
        presp
    DESTINATION
        bin
    )
