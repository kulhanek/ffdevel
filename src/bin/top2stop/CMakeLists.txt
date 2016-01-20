# ==============================================================================
# FFDevel CMake File
# ==============================================================================

# program objects --------------------------------------------------------------
SET(TOP2STOP_SRC
        MMTypes.cpp
        Top2STop.cpp
        Top2STopOptions.cpp
        )

# final build ------------------------------------------------------------------
ADD_EXECUTABLE(top2stop ${TOP2STOP_SRC})

TARGET_LINK_LIBRARIES(top2stop
        ${ASL_LIB_NAME}
        ${SCIMAFIC_CLIB_NAME}
        ${HIPOLY_LIB_NAME}
        )

INSTALL(TARGETS
            top2stop
        RUNTIME DESTINATION
            bin
        )

