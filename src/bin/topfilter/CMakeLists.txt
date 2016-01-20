# ==============================================================================
# FFDevel CMake File
# ==============================================================================

# program objects --------------------------------------------------------------
SET(TOPFILTER_SRC
        TopFilter.cpp
        TopFilterOptions.cpp
        )

# final build ------------------------------------------------------------------
ADD_EXECUTABLE(topfilter ${TOPFILTER_SRC})

TARGET_LINK_LIBRARIES(topfilter
        ${ASL_LIB_NAME}
        ${HIPOLY_LIB_NAME}
        )

INSTALL(TARGETS
            topfilter
        RUNTIME DESTINATION
            bin
        )

