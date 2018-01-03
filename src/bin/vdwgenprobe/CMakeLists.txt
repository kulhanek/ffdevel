# ==============================================================================
# FFDevel CMake File
# ==============================================================================

# program objects --------------------------------------------------------------
SET(VDWGENPROBE_SRC
        VDWGenProbe.cpp
        VDWGenProbeOptions.cpp
        Sphere.cpp
        )

# final build ------------------------------------------------------------------
ADD_EXECUTABLE(vdwgenprobe ${VDWGENPROBE_SRC})

TARGET_LINK_LIBRARIES(vdwgenprobe
        ${OPEN_BABEL_LIB}
        ${SCIMAFIC_CLIB_NAME}
        ${HIPOLY_LIB_NAME}
        )

INSTALL(TARGETS
            vdwgenprobe
        RUNTIME DESTINATION
            bin
        )

