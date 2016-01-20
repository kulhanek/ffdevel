# ==============================================================================
# FFDevel CMake File
# ==============================================================================

# program objects --------------------------------------------------------------
SET(GENROTORS_SRC
        GenRotors.cpp
        GenRotorsOptions.cpp
        )

# final build ------------------------------------------------------------------
ADD_EXECUTABLE(genrotors ${GENROTORS_SRC})
ADD_DEPENDENCIES(genrotors ffdevel_shared)

TARGET_LINK_LIBRARIES(genrotors
        ${ASL_LIB_NAME}
        ${OPEN_BABEL_LIB}
        ${HIPOLY_LIB_NAME}
        )

INSTALL(TARGETS
            genrotors
        RUNTIME DESTINATION
            bin
        )

