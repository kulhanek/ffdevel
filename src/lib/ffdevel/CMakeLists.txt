# ==============================================================================
# FFDevel CMake File
# ==============================================================================

# objects in library -----------------------------------------------------------
SET(FFDEVEL_SRC
        sizes.f90
        constants.f90
        variables.f90
        ver.f90
        utils.f90
        timers.F90
        parallel_dat.f90
        parallel.f90

        forcefield/disp_dat.f90
        forcefield/disp.f90
        forcefield/disp_control.f90
        forcefield/topology_dat.f90
        forcefield/topology.f90
        forcefield/geometry_dat.f90
        forcefield/geometry.f90
        forcefield/geometry_utils.f90
        forcefield/energy.f90
        forcefield/gradient.f90
        forcefield/gradient_utils.f90
        forcefield/hessian.f90
        forcefield/hessian_utils.f90
        forcefield/jacobian.f90

        forcefield/nbmodes/LJ.f90
        forcefield/nbmodes/12-DISPBJ.f90
        forcefield/nbmodes/EXP-DISPBJ.f90
        forcefield/nbmodes/EXP-DISPTT.f90

        parameters/errors/errors.f90
        parameters/errors/errors_dat.f90
        parameters/errors/errors_control.f90
        parameters/errors/errors_utils.f90

        parameters/errors/bonds_dat.f90
        parameters/errors/bonds.f90
        parameters/errors/bonds_control.f90
        parameters/errors/angles_dat.f90
        parameters/errors/angles.f90
        parameters/errors/angles_control.f90
        parameters/errors/dihedrals_dat.f90
        parameters/errors/dihedrals.f90
        parameters/errors/dihedrals_control.f90
        parameters/errors/impropers_dat.f90
        parameters/errors/impropers.f90
        parameters/errors/impropers_control.f90
        parameters/errors/nbdists_dat.f90
        parameters/errors/nbdists.f90
        parameters/errors/nbdists_control.f90
        parameters/errors/rmsd_dat.f90
        parameters/errors/rmsd.f90
        parameters/errors/rmsd_control.f90
        parameters/errors/ihess_dat.f90
        parameters/errors/ihess.f90
        parameters/errors/ihess_control.f90
        parameters/errors/energy_dat.f90
        parameters/errors/energy.f90
        parameters/errors/energy_control.f90
        parameters/errors/sapt_dat.f90
        parameters/errors/sapt.f90
        parameters/errors/sapt_control.f90
        parameters/errors/chrgpnl_dat.f90
        parameters/errors/chrgpnl.f90
        parameters/errors/chrgpnl_control.f90
        parameters/errors/zerograd_dat.f90
        parameters/errors/zerograd.f90
        parameters/errors/zerograd_control.f90
        parameters/errors/probe_dat.f90
        parameters/errors/probe.f90
        parameters/errors/probe_control.f90
        parameters/errors/pbpnl_dat.f90
        parameters/errors/pbpnl.f90
        parameters/errors/pbpnl_control.f90

        parameters/densoverlap_dat.f90
        parameters/densoverlap_db.f90
        parameters/densoverlap.f90
        parameters/densoverlap_control.f90
        parameters/buried_dat.f90
        parameters/buried.f90
        parameters/buried_control.f90
        parameters/xdm_dat.f90
        parameters/xdm.f90
        parameters/targetset_dat.f90
        parameters/targetset.f90
        parameters/targetset_control.f90
        parameters/parameters_dat.f90
        parameters/parameters.f90
        parameters/parameters_control.f90
        parameters/nb2nb.f90
        parameters/nb2nb_dat.f90
        parameters/nb2nb_control.f90

        optimizers/geoopt_dat.f90
        optimizers/geoopt_control.f90
        optimizers/geoopt.f90

        optimizers/ffopt_dat.f90
        optimizers/ffopt_control.f90
        optimizers/ffopt.f90
        optimizers/shark-ml1.cpp
        optimizers/shark-ml2.cpp

        dftd3/sizes.f90
        dftd3/api.f90
        dftd3/common.f90
        dftd3/core.f90
        dftd3/pars.f90
        dftd3/sizes.f90
        dftd3/mmd3.f90
        dftd3/mmd3_dat.f90
        )

# create static library --------------------------------------------------------
IF(LIBS_STATIC)
    ADD_LIBRARY(ffdevel_static STATIC ${FFDEVEL_SRC})
    SET_TARGET_PROPERTIES(ffdevel_static PROPERTIES
                            OUTPUT_NAME ffdevel
                            CLEAN_DIRECT_OUTPUT 1)
ENDIF(LIBS_STATIC)

# create shared library --------------------------------------------------------
IF(LIBS_SHARED)
    ADD_LIBRARY(ffdevel_shared SHARED ${FFDEVEL_SRC})
    SET_TARGET_PROPERTIES(ffdevel_shared PROPERTIES
                            OUTPUT_NAME ffdevel
                            CLEAN_DIRECT_OUTPUT 1
                            VERSION ${FFDEVEL_VERSION})
    TARGET_LINK_LIBRARIES(ffdevel_shared
            ${SHARK_LIB_NAME}
            ${FFTW_LIB_NAME}
            ${NLOPT_LIB_NAME}
            ${LBFGS_LIB_NAME}
            ${PRMFILE_FLIB_NAME}
            ${SCIMAFIC_FLIB_NAME}
            ${LAPACK_LIB_NAME}
            ${BLAS_LIB_NAME}
            )
    INSTALL(TARGETS
                ffdevel_shared
            LIBRARY DESTINATION
                lib
            )
ENDIF(LIBS_SHARED)
