# ==============================================================================
# FFDevel CMake File
# ==============================================================================

# program objects --------------------------------------------------------------
SET(TOPREPORT_SRC
        MMTypes.cpp
        TopReport.cpp
        TopReportOptions.cpp
        )

# final build ------------------------------------------------------------------
ADD_EXECUTABLE(topreport ${TOPREPORT_SRC})

TARGET_LINK_LIBRARIES(topreport
        ${ASL_LIB_NAME}
        ${HIPOLY_LIB_NAME}
        )

INSTALL(TARGETS
            topreport
        RUNTIME DESTINATION
            bin
        )

