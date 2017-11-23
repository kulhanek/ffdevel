# ==============================================================================
# FFDevel CMake File
# ==============================================================================

# program objects --------------------------------------------------------------
SET(VDWGENPROBE_SRC
        VDWGenProbe.cpp
        VDWGenProbeOptions.cpp
        )

# final build ------------------------------------------------------------------
ADD_EXECUTABLE(vdwgenprobe ${VDWGENPROBE_SRC})

TARGET_LINK_LIBRARIES(vdwgenprobe
        ${SCIMAFIC_CLIB_NAME}
        ${HIPOLY_LIB_NAME}
        )

INSTALL(TARGETS
            vdwgenprobe
        RUNTIME DESTINATION
            bin
        )

