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

program ffdev_params_program

    use ffdev_sizes
    use ffdev_utils
    use ffdev_constants
    use prmfile
    use ffdev_targetset_control
    use ffdev_parameters
    use ffdev_parameters_control

    implicit none
    character(len=MAX_PATH) :: ctrlname      ! input control file name
    type(PRMFILE_TYPE)      :: fin
    logical                 :: rst
    ! --------------------------------------------------------------------------

    call ffdev_utils_header('FF Parameters')

    ! test number of input arguments
    if( command_argument_count() .ne. 1 ) then
        call print_usage()
        call ffdev_utils_exit(DEV_OUT,1,'Incorrect number of arguments was specified (one expected)!')
    end if

    call get_command_argument(1, ctrlname)

    ! process control file -----------------------------------------------------
    write(DEV_OUT,*)
    call ffdev_utils_heading(DEV_OUT,'Reading control file', ':')
    write(DEV_OUT,'(a,a)') 'Control file name : ',trim(ctrlname)

    call prmfile_init(fin)

    if( .not. prmfile_read(fin,ctrlname) ) then
        call ffdev_utils_exit(DEV_OUT,1,'Specified control file cannot be opened!')
    end if

    ! read sections
    call ffdev_targetset_ctrl(fin,.true.)

    ! generate parameters from target topologies
    call ffdev_parameters_init()

    ! check if everything was read
    if( prmfile_count_ulines(fin,'TARGETS') .ne. 0 ) then
        write(DEV_OUT,*)
        call prmfile_dump_group(fin,DEV_OUT,'TARGETS',.true.)
        call ffdev_utils_exit(DEV_OUT,1,'Unprocessed lines found in the control file!')
    end if

    ! do we have identities in the {MAIN} group?
    rst = prmfile_open_group(fin,'MAIN')
    if( rst ) then
        call ffdev_parameters_ctrl_identities(fin)
        call ffdev_parameters_print_parameters()
    end if

    ! release the file
    call prmfile_clear(fin)

    call ffdev_utils_footer('FF Parameters')

100 format('Control file : ',A)

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

10 format('    ffparams <controlfile>')

end subroutine print_usage

!===============================================================================

end program ffdev_params_program
