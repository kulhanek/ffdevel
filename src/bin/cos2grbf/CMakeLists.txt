# ==============================================================================
# FFDevel CMake File
# ==============================================================================

# program objects --------------------------------------------------------------
SET(COS2GRBF_SRC
        ../top2stop/MMTypes.cpp
        Cos2GRBF.cpp
        Cos2GRBFOptions.cpp
        )

# final build ------------------------------------------------------------------
ADD_EXECUTABLE(cos2grbf ${COS2GRBF_SRC})

TARGET_LINK_LIBRARIES(cos2grbf
        ${ASL_LIB_NAME}
        ${SCIMAFIC_CLIB_NAME}
        ${HIPOLY_LIB_NAME}
        )

INSTALL(TARGETS
            cos2grbf
        RUNTIME DESTINATION
            bin
        )

