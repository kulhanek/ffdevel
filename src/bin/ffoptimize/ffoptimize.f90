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
    use prmfile
    use ffdev_targetset_control
    use ffdev_parameters
    use ffdev_parameters_dat
    use ffdev_parameters_control
    use ffdev_targetset
    use ffdev_mmd3
    use ffdev_mmd3_dontrol

    implicit none
    character(len=MAX_PATH)     :: ctrlname      ! input control file name
    type(PRMFILE_TYPE)          :: fin
    logical                     :: rst
    character(PRMFILE_MAX_PATH) :: string
    integer                     :: i
    ! --------------------------------------------------------------------------

    call ffdev_utils_header('FF Optimize')

    ! test number of input arguments
    if( command_argument_count() .ne. 1 ) then
        call print_usage()
        call ffdev_utils_exit(DEV_OUT,1,'Incorrect number of arguments was specified (one expected)!')
    end if

    call get_command_argument(1, ctrlname)
    
    ! open dev null
    call ffdev_utils_open(DEV_NULL,'/dev/null','O')

    ! process control file -----------------------------------------------------
    write(DEV_OUT,*)
    call ffdev_utils_heading(DEV_OUT,'Reading control file', ':')
    write(DEV_OUT,'(a,a)') 'Control file name : ',trim(ctrlname)

    call prmfile_init(fin)

    if( .not. prmfile_read(fin,ctrlname) ) then
        call ffdev_utils_exit(DEV_OUT,1,'Specified control file cannot be opened!')
    end if

    ! read files
    if( prmfile_open_group(fin,'MAIN') ) then
        call ffdev_parameters_ctrl_control(fin)
        call ffdev_parameters_ctrl_files(fin)
    end if

    ! init mmd3
    call ffdev_mmd3_init
    call ffdev_mmd3_ctrl(fin)

    ! read sections
    call ffdev_targetset_ctrl(fin,.false.)

    ! generate parameters from target topologies
    call ffdev_parameters_init()

    ! check if everything was read for TARGETS
    if( prmfile_count_ulines(fin,'TARGETS') .ne. 0 ) then
        write(DEV_OUT,*)
        call prmfile_dump_group(fin,DEV_OUT,'TARGETS',.true.)
        call ffdev_utils_exit(DEV_OUT,1,'Unprocessed lines found in the control file!')
    end if

    ! finalize topologies in sets
    call ffdev_targetset_init_pts()

    ! do fake input file processing
    rst = prmfile_first_group(fin)
    i = 1
    do while( rst )

        ! open set section
        if( .not. prmfile_get_group_name(fin,string) ) then
            call ffdev_utils_exit(DEV_OUT,1,'Unable to get group name!')
        end if

        if( string .eq. 'PROGRAM' ) then
            write(DEV_OUT,*)
            write(string,110) i
            call ffdev_utils_heading(DEV_OUT,trim(string), ':')

            call execute_program_fake(fin)

            i = i + 1
        end if

        if( string .eq. 'FFERROR' ) then
            call execute_fferror(fin,.true.)
        end if

        if( string .eq. 'FFMANIP' ) then
            call execute_ffmanip(fin,.true.)
        end if

        rst = prmfile_next_group(fin)
    end do

    ! check if everything was read
    if( prmfile_count_ulines(fin) .ne. 0 ) then
        write(DEV_OUT,*)
        call prmfile_dump(fin,DEV_OUT,.true.)
        call ffdev_utils_exit(DEV_OUT,1,'Unprocessed lines found in the control file!')
    end if

    ! start optimization programs ----------------
    write(DEV_OUT,*)
    call ffdev_utils_heading(DEV_OUT,'==========================', '!')
    call ffdev_utils_heading(DEV_OUT,'Starting real optimization', '!')
    call ffdev_utils_heading(DEV_OUT,'==========================', '!')

    ! process optimization programs
    rst = prmfile_first_group(fin)
    i = 1
    do while( rst )

        ! open set section
        if( .not. prmfile_get_group_name(fin,string) ) then
            call ffdev_utils_exit(DEV_OUT,1,'Unable to get group name!')
        end if

    ! fferror -------------------------------------
        if( string .eq. 'FFERROR' ) then
            call execute_fferror(fin,.false.)
        end if

    ! manipff -------------------------------------
        if( string .eq. 'FFMANIP' ) then
            call execute_ffmanip(fin,.false.)
        end if

    ! program -------------------------------------
        if( string .eq. 'PROGRAM' ) then
            write(DEV_OUT,*)
            write(string,110) i
            call ffdev_utils_heading(DEV_OUT,trim(string), ':')

            call execute_program(fin)

            i = i + 1
        end if

        rst = prmfile_next_group(fin)
    end do

    ! check if everything was read
    if( prmfile_count_ulines(fin) .ne. 0 ) then
        write(DEV_OUT,*)
        call prmfile_dump(fin,DEV_OUT,.true.)
        call ffdev_utils_exit(DEV_OUT,1,'Unprocessed lines found in the control file!')
    end if

    ! release the file
    call prmfile_clear(fin)

    ! finalize
    write(DEV_OUT,*)
    call ffdev_utils_heading(DEV_OUT,'Finalize', ':')

    ! be sure that optimized parameters are in all topologies
    ! (which is not case if numerical gradient of error function is employed)
    call ffdev_parameters_to_tops()

    ! save topologies if requested
    call ffdev_targetset_save_final_stops()

    ! save geometries if requested
    call ffdev_targetset_save_final_pts()

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

    call ffdev_utils_footer('FF Optimize')

100 format('Control file : ',A)
110 format('FF optimization program #',I2.2)
120 format('Final ffdevel parameters        = ',A)
130 format('Final AMBER parameter file      = ',A)

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
! subroutine:  execute_program
!===============================================================================

subroutine execute_program(grpin)

    use ffdev_parameters_control
    use ffdev_ffopt_control
    use ffdev_ffopt

    implicit none
    type(PRMFILE_TYPE)  :: grpin
    ! --------------------------------------------------------------------------

    ! load setup
    call ffdev_parameters_ctrl_identities(grpin)
    call ffdev_parameters_ctrl_realms(grpin)
    call ffdev_ffopt_ctrl_minimize(grpin)

    ! print final parameter list
    call ffdev_parameters_print_parameters()

    ! optimize
    call ffdev_ffopt_run

    ! final results
    call ffdev_parameters_print_parameters()

    return

end subroutine execute_program

!===============================================================================
! subroutine:  execute_program_fake
!===============================================================================

subroutine execute_program_fake(grpin)

    use ffdev_parameters_control
    use ffdev_ffopt_control
    use ffdev_ffopt

    implicit none
    type(PRMFILE_TYPE)  :: grpin
    ! --------------------------------------------------------------------------

    ! load setup
    call ffdev_parameters_ctrl_identities(grpin)
    call ffdev_parameters_ctrl_realms(grpin)
    call ffdev_ffopt_ctrl_minimize(grpin)

    return

end subroutine execute_program_fake

!===============================================================================
! subroutine:  execute_ffmanip
!===============================================================================

subroutine execute_ffmanip(grpin,noexec)

    use ffdev_parameters_control
    use ffdev_ffopt_control

    implicit none
    type(PRMFILE_TYPE)  :: grpin
    logical             :: noexec
    ! --------------------------------------------------------------------------

    if( .not. noexec ) then
        call ffdev_parameters_print_parameters()
    end if

    ! load and execute setup
    call ffdev_parameters_ctrl_ffmanip(grpin,noexec)

    if( .not. noexec ) then
        call ffdev_parameters_print_parameters()
    end if

    return

end subroutine execute_ffmanip

!===============================================================================
! subroutine:  execute_fferror
!===============================================================================

subroutine execute_fferror(grpin,noexec)

    use ffdev_errors_control

    implicit none
    type(PRMFILE_TYPE)  :: grpin
    logical             :: noexec
    ! --------------------------------------------------------------------------

    ! load and execute error setup
    call ffdev_errors_ctrl(grpin)

    return

end subroutine execute_fferror

!===============================================================================

end program ffdev_optimize_program
