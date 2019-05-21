# ==============================================================================
# FFDevel CMake File
# ==============================================================================

INCLUDE_DIRECTORIES(lib/ffdevel SYSTEM)
INCLUDE_DIRECTORIES(lib/ffdevel/optimizers SYSTEM)
INCLUDE_DIRECTORIES(lib/numgrid SYSTEM)

ADD_SUBDIRECTORY(lib)
ADD_SUBDIRECTORY(bin)



