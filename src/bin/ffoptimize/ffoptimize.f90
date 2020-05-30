! ==============================================================================
! This file is part of FFDevel.
!    Copyright (C) 2013 Petr Kulhanek, kulhanek@chemi.muni.cz
!
! FFDevel is free software: you can redistribute it and/or modify it under
! the terms of the GNU General Public License as published by the Free Software
! Foundation, either version 3 of the License, or (at your option) any later
! version.
!
! FFDevel is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
! FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with FFDevel. If not, see <http://www.gnu.org/licenses/>.
! ==============================================================================

program ffdev_optimize_program

    use ffdev_sizes
    use ffdev_utils
    use ffdev_constants
    use ffdev_variables
    use prmfile
    use ffdev_targetset_control
    use ffdev_parameters
    use ffdev_parameters_dat
    use ffdev_parameters_control
    use ffdev_targetset
    use ffdev_errors
    use ffdev_errors_dat
    use ffdev_xdm
    use ffdev_mmd3
    use ffdev_timers
!$ use omp_lib

    implicit none
    character(len=MAX_PATH)     :: ctrlname     ! input control file name
    type(PRMFILE_TYPE)          :: fin
    type(PRMFILE_TYPE)          :: tmpfin
    logical                     :: rst
    character(PRMFILE_MAX_PATH) :: string
    integer                     :: i, ncpu
    ! --------------------------------------------------------------------------

    call ffdev_utils_header('FF Optimize')

    call ffdev_timers_init_top
    call ffdev_timers_init

    call ffdev_timers_start_timer(FFDEV_INITIALIZATION_TIMER)

    ! test number of input arguments
    if( command_argument_count() .ne. 1 ) then
        call print_usage()
        call ffdev_utils_exit(DEV_ERR,1,'Incorrect number of arguments was specified (one expected)!')
    end if

    call get_command_argument(1, ctrlname)

    ! open dev null
    call ffdev_utils_open(DEV_NULL,'/dev/null','O')

! ==============================================================================
! INITIAL PART - Initial setup + TARGETSETS
! ==============================================================================

    ! default setup for subsystems ---------------------------------------------
    call prmfile_init(tmpfin)
    write(DEV_OUT,*)
    call ffdev_utils_heading(DEV_OUT,'Default setup of subsystems', ':')
    call ffdev_parameters_ctrl_nbsetup(tmpfin,.false.)
    call execute_mmopt(tmpfin,.false.)

    ! process control file -----------------------------------------------------
    write(DEV_OUT,*)
    call ffdev_utils_heading(DEV_OUT,'Reading control file', ':')
    write(DEV_OUT,'(a,a)') 'Control file name : ',trim(ctrlname)

    call prmfile_init(fin)

    if( .not. prmfile_read(fin,ctrlname) ) then
        call ffdev_utils_exit(DEV_ERR,1,'Specified control file cannot be opened!')
    end if

    ! read configuration
    if( prmfile_open_group(fin,'FFPARAMS') ) then
        write(DEV_OUT,*)
        call ffdev_utils_heading(DEV_OUT,'{FFPARAMS}',':')
        call ffdev_parameters_ctrl_control(fin)
        call ffdev_parameters_ctrl_files(fin)
        call ffdev_parameters_ctrl_grbf2cos(fin)
    end if

    ! read target sets
    call ffdev_targetset_ctrl(fin,.false.)

    ! check if everything was read for TARGETS
    if( prmfile_count_ulines(fin,'TARGETS') .ne. 0 ) then
        write(DEV_OUT,*)
        call prmfile_dump_group(fin,DEV_OUT,'TARGETS',.true.)
        call ffdev_utils_exit(DEV_ERR,1,'Unprocessed lines found in the control file!')
    end if

    ! generate parameters from target topologies
    call ffdev_parameters_init()

    ! finalize topologies in sets
    call ffdev_targetset_init_pts()

    ! run mmd3
    call ffdev_mmd3_init
    call ffdev_mmd3_run_stat()

    ! run XDM stat if data available
    call ffdev_xdm_run_stat()

! ==============================================================================
! INITIAL PART - check for syntax errors in control file
! ==============================================================================

    ! start optimization programs ----------------
    write(DEV_OUT,*)
    call ffdev_utils_heading(DEV_OUT,'================================', '!')
    call ffdev_utils_heading(DEV_OUT,'Checking control file for errors', '!')
    call ffdev_utils_heading(DEV_OUT,'================================', '!')

    if( Verbosity .ge. DEV_VERBOSITY_FULL ) then
        DEV_OUT = DEV_STD_OUTPUT
    else
        DEV_OUT = DEV_NULL
    end if

    ! reset initial setup
    call ffdev_targetset_ctrl_optgeo_set_default()
    call ffdev_parameters_disable_all_realms()
    call ffdev_errors_init_all()

    ! do fake input file processing
    rst = prmfile_first_group(fin)
    i = 1
    do while( rst )

        ! open set section
        if( .not. prmfile_get_group_name(fin,string) ) then
            call ffdev_utils_exit(DEV_ERR,1,'Unable to get group name!')
        end if

        if( string .eq. 'CHECKGRD' ) then
            call execute_check_gradient(fin,.false.)
        end if

        if( string .eq. 'FFERROR' ) then
            call execute_fferror(fin,.false.)
        end if

        if( string .eq. 'FFMANIP' ) then
            call execute_ffmanip(fin,.false.)
        end if

        if( string .eq. 'MMOPT' ) then
            call execute_mmopt(fin,.false.)
        end if

        if( string .eq. 'FFOPT' ) then
            write(DEV_OUT,*)
            write(string,110) i
            call ffdev_utils_heading(DEV_OUT,trim(string), ':')

            call execute_ffopt(fin,i,.false.)

            i = i + 1
        end if

        if( string .eq. 'FFEVAL' ) then
            write(DEV_OUT,*)
            write(string,110) i
            call ffdev_utils_heading(DEV_OUT,trim(string), ':')

            call execute_ffeval(fin,i,.false.)

            i = i + 1
        end if

        rst = prmfile_next_group(fin)
    end do

    call ffdev_timers_stop_timer(FFDEV_INITIALIZATION_TIMER)

    ! restore output channel
    DEV_OUT = DEV_STD_OUTPUT

    ! check if everything was read
    if( prmfile_count_ulines(fin) .ne. 0 ) then
        write(DEV_ERR,*)
        call prmfile_dump(fin,DEV_ERR,.true.)
        call ffdev_utils_exit(DEV_ERR,1,'Unprocessed lines found in the control file!')
    else
        write(DEV_OUT,*) '>>> INFO: The control file seems to be OK!'
    end if

    ! start optimization programs ----------------
    write(DEV_OUT,*)
    call ffdev_utils_heading(DEV_OUT,'==========================', '!')
    call ffdev_utils_heading(DEV_OUT,'Starting real optimization', '!')
    call ffdev_utils_heading(DEV_OUT,'==========================', '!')

    ncpu = 1
    !$ ncpu = omp_get_max_threads()
    write(*,'(A,I3)') '>>> INFO: Number of threads = ',ncpu

    ! reset initial setup
    call ffdev_targetset_ctrl_optgeo_set_default()
    call ffdev_parameters_disable_all_realms()
    call ffdev_errors_init_all()

    ! initialize
    write(DEV_OUT,*)
    call ffdev_utils_heading(DEV_OUT,'Initialize', ':')

    ! FIXME this must be done in FFOPT after all is setup, for exemple in FFMANIP
!    ! calculate initial data
!    errors_calc_ene = .true. ! we need at least energies for ffdev_targetset_save_initial_drvs
!    call ffdev_targetset_calc_all()

!    ! save initial driving data if requested
!    ! FIXME - the geometry cannot be never optimized as default setup
!    ! set in ffdev_targetset_ctrl_optgeo_set_default, see above?
!    ! the only solution is to enable geometry optimization directly for the given [set]
!    call ffdev_targetset_save_initial_drvs()

    ! process optimization programs
    rst = prmfile_first_group(fin)
    i = 1
    do while( rst )

        ! open set section
        if( .not. prmfile_get_group_name(fin,string) ) then
            call ffdev_utils_exit(DEV_ERR,1,'Unable to get group name!')
        end if

    ! check gradient ------------------------------
        if( string .eq. 'CHECKGRD' ) then
            call execute_check_gradient(fin,.true.)
        end if

    ! fferror -------------------------------------
        if( string .eq. 'FFERROR' ) then
            call execute_fferror(fin,.true.)
        end if

    ! ffmanip -------------------------------------
        if( string .eq. 'FFMANIP' ) then
            call execute_ffmanip(fin,.true.)
        end if

    ! mmopt ---------------------------------------
        if( string .eq. 'MMOPT' ) then
            call execute_mmopt(fin,.true.)
        end if

    ! ffopt ---------------------------------------
        if( string .eq. 'FFOPT' ) then
            write(DEV_OUT,*)
            write(string,110) i
            call ffdev_utils_heading(DEV_OUT,trim(string), ':')

            call execute_ffopt(fin,i,.true.)

            i = i + 1
        end if

    ! ffeval - single point error evaluation
        if( string .eq. 'FFEVAL' ) then
            write(DEV_OUT,*)
            write(string,110) i
            call ffdev_utils_heading(DEV_OUT,trim(string), ':')

            call execute_ffeval(fin,i,.true.)

            i = i + 1
        end if

        rst = prmfile_next_group(fin)
    end do

    ! check if everything was read
    if( prmfile_count_ulines(fin) .ne. 0 ) then
        write(DEV_OUT,*)
        call prmfile_dump(fin,DEV_OUT,.true.)
        call ffdev_utils_exit(DEV_ERR,1,'Unprocessed lines found in the control file!')
    end if

    ! release the file
    call prmfile_clear(fin)

    ! finalize
    write(DEV_OUT,*)
    call ffdev_utils_heading(DEV_OUT,'Finalize', ':')

    ! be sure that optimized parameters are in all topologies
    ! (which is not case if numerical gradient of error function is employed)
    call ffdev_parameters_to_tops()

    ! end recalculate sets
    call ffdev_targetset_calc_all()

    ! save topologies if requested
    call ffdev_targetset_save_final_stops()

    ! save drivings if requested
    call ffdev_targetset_save_final_drvs()

    ! save geometries if requested
    call ffdev_targetset_save_final_pts()
    call ffdev_targetset_save_final_xyzr()

    ! save optimized parameters
    if( len_trim(OutParamFileName) .ne. 0 ) then
        write(DEV_OUT,120) trim(OutParamFileName)
        call ffdev_parameters_save(OutParamFileName)
    end if

    ! save final amber parameter file
    if( len_trim(OutAmberPrmsFileName) .ne. 0 ) then
        write(DEV_OUT,120) trim(OutAmberPrmsFileName)
        call ffdev_parameters_save_amber(OutAmberPrmsFileName)
    end if

    call ffdev_timers_finalize(.true.)

    call ffdev_utils_footer('FF Optimize')

110 format('FF optimization program #',I3.3)
120 format('Final ffdevel parameters        = ',A)

contains

!===============================================================================
! subroutine:  print_usage
!===============================================================================

subroutine print_usage()

    implicit none
    ! --------------------------------------------------------------------------

    write(DEV_OUT,*)
    write(DEV_OUT,'(/,a,/)') '=== [usage] ===================================================================='
    write(DEV_OUT,10)
    write(DEV_OUT,*)

    return

10 format('    ffoptimize <controlfile>')

end subroutine print_usage

!===============================================================================
! subroutine:  execute_ffopt
!===============================================================================

subroutine execute_ffopt(grpin,pid,exec)

    use ffdev_parameters_control
    use ffdev_ffopt_control
    use ffdev_ffopt

    implicit none
    type(PRMFILE_TYPE)  :: grpin
    integer             :: pid
    logical             :: exec
    ! --------------------------------------------------------------------------

    CurrentProgID = pid

    ! load setup
    call ffdev_parameters_ctrl_identities(grpin)
    call ffdev_parameters_ctrl_realms(grpin)
    call ffdev_targetset_ctrl_optgeo(grpin)
    call ffdev_ffopt_ctrl_minimize(grpin)

    if( exec ) then
        ! print final parameter list
        call ffdev_parameters_print_parameters(PARAMS_SUMMARY_INITIAL)

        ! optimize
        call ffdev_ffopt_run

        ! final results
        call ffdev_parameters_print_parameters(PARAMS_SUMMARY_OPTIMIZED)
    end if

    return

end subroutine execute_ffopt

!===============================================================================
! subroutine:  execute_ffeval
!===============================================================================

subroutine execute_ffeval(grpin,pid,exec)

    use ffdev_parameters_control
    use ffdev_ffopt_control
    use ffdev_ffopt

    implicit none
    type(PRMFILE_TYPE)  :: grpin
    integer             :: pid
    logical             :: exec
    ! --------------------------------------------------------------------------

    CurrentProgID = pid

    ! load setup
    call ffdev_parameters_ctrl_identities(grpin)
    call ffdev_parameters_ctrl_realms(grpin)
    call ffdev_targetset_ctrl_optgeo(grpin)

    if( exec ) then
        ! print final parameter list
        call ffdev_parameters_print_parameters(PARAMS_SUMMARY_FULL)

        ! get single point error
        call ffdev_ffopt_single_point
    end if

    return

end subroutine execute_ffeval

!===============================================================================
! subroutine:  execute_mmopt
!===============================================================================

subroutine execute_mmopt(grpin,exec)

    use ffdev_parameters_control
    use ffdev_geoopt_control

    implicit none
    type(PRMFILE_TYPE)  :: grpin
    logical             :: exec
    ! --------------------------------------------------------------------------

    write(DEV_OUT,*)
    call ffdev_utils_heading(DEV_OUT,'MMOPT', ':')

    ! load setup
    call ffdev_geoopt_ctrl_minimize(grpin)

    return

end subroutine execute_mmopt

!===============================================================================
! subroutine:  execute_ffmanip
!===============================================================================

subroutine execute_ffmanip(grpin,exec)

    use ffdev_parameters_control
    use ffdev_ffopt_control

    implicit none
    type(PRMFILE_TYPE)  :: grpin
    logical             :: exec
    ! --------------------------------------------------------------------------

    if( exec ) then
        call ffdev_parameters_print_parameters(PARAMS_SUMMARY_FULL)
    end if

    ! load and execute setup
    call ffdev_parameters_ctrl_ffmanip(grpin,exec)

    if( exec ) then
        call ffdev_parameters_print_parameters(PARAMS_SUMMARY_MODIFIED)
    end if

    return

end subroutine execute_ffmanip

!===============================================================================
! subroutine:  execute_fferror
!===============================================================================

subroutine execute_fferror(grpin,exec)

    use ffdev_errors_control

    implicit none
    type(PRMFILE_TYPE)  :: grpin
    logical             :: exec
    ! --------------------------------------------------------------------------

    ! load and execute error setup
    call ffdev_errors_ctrl(grpin)

    return

end subroutine execute_fferror

!===============================================================================
! subroutine:  execute_check_gradient
!===============================================================================

subroutine execute_check_gradient(grpin,exec)

    use ffdev_targetset_dat
    use ffdev_gradient
    use ffdev_gradient_utils
    use ffdev_geometry
    use ffdev_geometry_utils

    implicit none
    type(PRMFILE_TYPE)  :: grpin
    logical             :: exec
    type(GEOMETRY)      :: ageo,ngeo
    ! --------------------------------------------
    integer             :: i,j
    ! --------------------------------------------------------------------------

    if( .not. exec ) return

    do i=1,nsets
        do j=1,sets(i)%ngeos
            write(DEV_OUT,*)
            write(DEV_OUT,1) i,j
            write(DEV_OUT,'(A)') 'Analytical gradient ...'
            call ffdev_geometry_init(ageo)
            call ffdev_geometry_copy(ageo,sets(i)%geo(j))
            call ffdev_gradient_allocate(ageo)
            call ffdev_gradient_all(sets(i)%top,ageo)
            call ffdev_geometry_info_ene(ageo)

            write(DEV_OUT,'(A)') 'Numerical gradient ...'
            call ffdev_geometry_init(ngeo)
            call ffdev_geometry_copy(ngeo,sets(i)%geo(j))
            call ffdev_gradient_allocate(ngeo)
            call ffdev_gradient_num_all(sets(i)%top,ngeo)
            call ffdev_geometry_info_ene(ngeo)

            if( abs(ageo%total_ene - ngeo%total_ene) .gt. 1e-6 ) then
                call ffdev_utils_exit(DEV_ERR,1,'Analytical and numerical total energies do not match!')
            end if

            if( .not. ffdev_gradient_test(ageo,ngeo,1.0d-3) ) then
                write(DEV_OUT,*)
                call ffdev_utils_heading(DEV_OUT,'Analytical FF Gradient','=')
                call ffdev_gradient_print(DEV_OUT,sets(i)%top,ageo)
                    write(DEV_OUT,*)
                call ffdev_utils_heading(DEV_OUT,'Numerical FF Gradient','=')
                call ffdev_gradient_print(DEV_OUT,sets(i)%top,ngeo)

                write(DEV_OUT,*)
                call ffdev_utils_heading(DEV_OUT,'Difference in Gradients','=')
                ngeo%grd = ngeo%grd - ageo%grd
                call ffdev_gradient_print(DEV_OUT,sets(i)%top,ngeo)

                call ffdev_utils_exit(DEV_ERR,1,'Analytical and numerical gradients do not match!')
            else
                write(DEV_OUT,'(A)') 'Analytical and numerical gradients match each other ...'
            end if
        end do
    end do

  1 format('=== [SET ',I3.3,']/[GEO ',I3.3,'] ========================================================')

end subroutine execute_check_gradient

!===============================================================================

end program ffdev_optimize_program
