# ==============================================================================
# FFDevel CMake File
# ==============================================================================

CMAKE_MINIMUM_REQUIRED(VERSION 2.8 FATAL_ERROR)
#SET(CMAKE_VERBOSE_MAKEFILE ON)

# define project name ----------------------------------------------------------
PROJECT(FFDevel)
SET(FFDEVEL_VERSION 1.0.x)

# what languages we use --------------------------------------------------------
ENABLE_LANGUAGE(CXX)
ENABLE_LANGUAGE(Fortran)

# ==============================================================================
# conditional setup ------------------------------------------------------------
# ==============================================================================

SET(LIBS_SHARED     ON  CACHE BOOL "Should the static version of xdyn libraries be built?")
SET(LIBS_STATIC     OFF CACHE BOOL "Should the static version of prmfile library be built?")

# ==============================================================================
# external dependencies --------------------------------------------------------
# ==============================================================================

IF(NOT DEFINED COMPOSITE_PROJECT)
    SET(DEVELOPMENT_ROOT "$ENV{DEVELOPMENT_ROOT}")
ENDIF(NOT DEFINED COMPOSITE_PROJECT)

IF(NOT DEVELOPMENT_ROOT)
    MESSAGE(FATAL_ERROR "DEVELOPMENT_ROOT variable is not set!")
ENDIF(NOT DEVELOPMENT_ROOT)

# ==============================================================================
# external libraries -----------------------------------------------------------
# ==============================================================================

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${DEVELOPMENT_ROOT}/projects/cmake-modules/1.0/share/modules")

# boost =======================
FIND_PACKAGE(MyBoost REQUIRED)
INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS} SYSTEM)

# HIPOLY ====================
SET(HIPOLY_ROOT ${DEVELOPMENT_ROOT}/projects/hipoly/1.0)
INCLUDE_DIRECTORIES(${HIPOLY_ROOT}/src/lib/hipoly SYSTEM)
INCLUDE_DIRECTORIES(${HIPOLY_ROOT}/src/lib/hipoly/compress SYSTEM)
INCLUDE_DIRECTORIES(${HIPOLY_ROOT}/src/lib/hipoly/error SYSTEM)
INCLUDE_DIRECTORIES(${HIPOLY_ROOT}/src/lib/hipoly/event SYSTEM)
INCLUDE_DIRECTORIES(${HIPOLY_ROOT}/src/lib/hipoly/file SYSTEM)
INCLUDE_DIRECTORIES(${HIPOLY_ROOT}/src/lib/hipoly/list SYSTEM)
INCLUDE_DIRECTORIES(${HIPOLY_ROOT}/src/lib/hipoly/option SYSTEM)
INCLUDE_DIRECTORIES(${HIPOLY_ROOT}/src/lib/hipoly/string SYSTEM)
INCLUDE_DIRECTORIES(${HIPOLY_ROOT}/src/lib/hipoly/stream SYSTEM)
INCLUDE_DIRECTORIES(${HIPOLY_ROOT}/src/lib/hipoly/thread SYSTEM)
INCLUDE_DIRECTORIES(${HIPOLY_ROOT}/src/lib/hipoly/time SYSTEM)
INCLUDE_DIRECTORIES(${HIPOLY_ROOT}/src/lib/hipoly/color SYSTEM)
INCLUDE_DIRECTORIES(${HIPOLY_ROOT}/src/lib/hipoly/xml SYSTEM)
LINK_DIRECTORIES(${HIPOLY_ROOT}/lib)
SET(HIPOLY_LIB_NAME hipoly)

# SCIMAFIC ====================
SET(SCIMAFIC_ROOT ${DEVELOPMENT_ROOT}/projects/scimafic/1.0)
INCLUDE_DIRECTORIES(${SCIMAFIC_ROOT}/src/lib/cscimafic SYSTEM)
INCLUDE_DIRECTORIES(${SCIMAFIC_ROOT}/src/lib/cscimafic/files SYSTEM)
INCLUDE_DIRECTORIES(${SCIMAFIC_ROOT}/src/lib/cscimafic/maths SYSTEM)
INCLUDE_DIRECTORIES(${SCIMAFIC_ROOT}/src/lib/cscimafic/chemistry SYSTEM)
INCLUDE_DIRECTORIES(${SCIMAFIC_ROOT}/src/lib/cscimafic/physics SYSTEM)
INCLUDE_DIRECTORIES(${SCIMAFIC_ROOT}/src/lib/fscimafic SYSTEM)
LINK_DIRECTORIES(${SCIMAFIC_ROOT}/lib)
SET(SCIMAFIC_CLIB_NAME cscimafic)
SET(SCIMAFIC_FLIB_NAME fscimafic)

# ASL =========================
SET(ASL_ROOT ${DEVELOPMENT_ROOT}/projects/asl/2.0)
INCLUDE_DIRECTORIES(${ASL_ROOT}/src/lib/asl SYSTEM)
INCLUDE_DIRECTORIES(${ASL_ROOT}/src/lib/asl/topology SYSTEM)
INCLUDE_DIRECTORIES(${ASL_ROOT}/src/lib/asl/restart SYSTEM)
INCLUDE_DIRECTORIES(${ASL_ROOT}/src/lib/asl/trajectory SYSTEM)
INCLUDE_DIRECTORIES(${ASL_ROOT}/src/lib/asl/mask SYSTEM)
LINK_DIRECTORIES(${ASL_ROOT}/lib)
SET(ASL_LIB_NAME asl)

# PRMFile ===================
SET(PRMFILE_ROOT ${DEVELOPMENT_ROOT}/projects/prmfile/1.0)
INCLUDE_DIRECTORIES(${PRMFILE_ROOT}/src/lib/fprmfile SYSTEM)
LINK_DIRECTORIES(${PRMFILE_ROOT}/lib)
SET(PRMFILE_FLIB_NAME fprmfile)

IF(NOT MKL_FOUND)
    # BLAS ======================
    SET(BLAS_ROOT ${DEVELOPMENT_ROOT}/thirdparty/blas/netlib)
    LINK_DIRECTORIES(${BLAS_ROOT}/lib)
    SET(BLAS_LIB_NAME fblasd)

    # LAPACK =====================
    SET(LAPACK_ROOT ${DEVELOPMENT_ROOT}/thirdparty/lapack/3.1.1)
    LINK_DIRECTORIES(${LAPACK_ROOT}/lib)
    SET(LAPACK_LIB_NAME flapackd)
ELSE(NOT MKL_FOUND)
    SET(LAPACK_LIB_NAME ${MKL_LIBS})
ENDIF(NOT MKL_FOUND)

# L-BFGS =====================
SET(LBFGS_ROOT ${DEVELOPMENT_ROOT}/thirdparty/lbfgs/1.0)
INCLUDE_DIRECTORIES(${LBFGS_ROOT}/src/lib/flbfgs)
LINK_DIRECTORIES(${LBFGS_ROOT}/lib)
SET(LBFGS_LIB_NAME flbfgs90)

# runlux =====================
SET(RUNLUX_ROOT ${DEVELOPMENT_ROOT}/thirdparty/ranlux/1.0)
LINK_DIRECTORIES(${RUNLUX_ROOT}/lib)
SET(RANLUX_LIB_NAME franlux)

# nlopt =====================
SET(NLOPT_ROOT ${DEVELOPMENT_ROOT}/thirdparty/nlopt/2.4.2)
INCLUDE_DIRECTORIES(${NLOPT_ROOT}/include)
LINK_DIRECTORIES(${NLOPT_ROOT}/lib)
SET(NLOPT_LIB_NAME nlopt)

# OpenBabel --------------------------------------
SET(OPEN_BABEL_ROOT ${DEVELOPMENT_ROOT}/thirdparty/obcore/2.4.0)
INCLUDE_DIRECTORIES(${OPEN_BABEL_ROOT}/include)
LINK_DIRECTORIES(${OPEN_BABEL_ROOT}/lib)
SET(OPEN_BABEL_LIB obcore)

# FFTW --------------------------------------
FIND_PACKAGE(FFTW3 REQUIRED)
INCLUDE_DIRECTORIES(${FFTW3_INCLUDE_DIRS} SYSTEM)
LINK_DIRECTORIES(${FFTW3_LIB_DIRS})
SET(FFTW_LIB_NAME fftw3)

# Shark-ML --------------------------------------
SET(SHARK_ROOT ${DEVELOPMENT_ROOT}/thirdparty/shark-ml/4.0.0)
INCLUDE_DIRECTORIES(${SHARK_ROOT}/include)
LINK_DIRECTORIES(${SHARK_ROOT}/lib)
SET(SHARK_LIB_NAME shark)

# ==============================================================================
# output directories -----------------------------------------------------------
# ==============================================================================

# define output directory for created libraries and module files ---------------
SET(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin)
SET(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR}/lib)
LINK_DIRECTORIES(${CMAKE_BINARY_DIR}/lib)

# ==============================================================================
# rpaths -----------------------------------------------------------------------
# ==============================================================================

IF(NOT DEFINED COMPOSITE_PROJECT)
    # use, i.e. don't skip the full RPATH for the build tree
    SET(CMAKE_SKIP_BUILD_RPATH  FALSE)

    # when building, don't use the install RPATH already
    # (but later on when installing)
    SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)

    # the RPATH to be used when installing
    SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

    # add the automatically determined parts of the RPATH
    # which point to directories outside the build tree to the install RPATH
    SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH FALSE)
ENDIF(NOT DEFINED COMPOSITE_PROJECT)

# ==============================================================================
# compiler setups --------------------------------------------------------------
# ==============================================================================

IF(NOT DEFINED COMPOSITE_PROJECT)
    # get the name of fortran compiler
    get_filename_component(FCOMP_NAME ${CMAKE_Fortran_COMPILER} NAME_WE)

    # default compiler setups - gfortran
    SET(CMAKE_Fortran_FLAGS_RELEASE "-O3 -march=core-avx2 -fopenmp")
    SET(CMAKE_Fortran_FLAGS_DEBUG   "-g -O0 -fcheck=all -ffpe-trap=invalid,zero,overflow -fopenmp -Wall")

    # ifort
    IF(FCOMP_NAME STREQUAL ifort)
        SET(CMAKE_Fortran_FLAGS_RELEASE "-O2")
        SET(CMAKE_Fortran_FLAGS_DEBUG   "-g -O0 -std03 -traceback -check all -warn all")
    ENDIF(FCOMP_NAME STREQUAL ifort)

    SET(CMAKE_BUILD_TYPE Release)
    #SET(CMAKE_BUILD_TYPE Debug)
ENDIF(NOT DEFINED COMPOSITE_PROJECT)

# ==============================================================================
# project subdirectories  ------------------------------------------------------
# ==============================================================================

ADD_SUBDIRECTORY(src)
