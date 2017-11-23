# ==============================================================================
# FFDevel CMake File
# ==============================================================================

# program objects --------------------------------------------------------------
SET(VDWPROBEFILTER_SRC
        VDWProbeFilter.cpp
        VDWProbeFilterOptions.cpp
        )

# final build ------------------------------------------------------------------
ADD_EXECUTABLE(vdwprobefilter ${VDWPROBEFILTER_SRC})

TARGET_LINK_LIBRARIES(vdwprobefilter
        ${SCIMAFIC_CLIB_NAME}
        ${HIPOLY_LIB_NAME}
        )

INSTALL(TARGETS
            vdwprobefilter
        RUNTIME DESTINATION
            bin
        )

