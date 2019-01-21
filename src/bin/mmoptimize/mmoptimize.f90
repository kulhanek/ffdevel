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

program ffdev_mmoptimize

    use ffdev_sizes
    use ffdev_utils
    use ffdev_constants
    use ffdev_topology
    use ffdev_geometry
    use prmfile
    use ffdev_geoopt_dat
    use ffdev_geoopt_control
    use ffdev_gradient_utils
    use ffdev_geoopt

    implicit none
    character(len=MAX_PATH) :: ctrlname      ! input control file name
    type(PRMFILE_TYPE)      :: fin
    type(TOPOLOGY)          :: top
    type(GEOMETRY)          :: geo
    ! --------------------------------------------------------------------------

    call ffdev_utils_header('MM Optimize')

    ! check if control file was provided
    if( (command_argument_count() .ne. 1) .and. (command_argument_count() .ne. 2) &
        .and. (command_argument_count() .ne. 3) .and. (command_argument_count() .ne. 4) ) then
        call print_usage
        call ffdev_utils_exit(DEV_OUT,1,'No input file specified on the command line!')
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

    if( .not. prmfile_open_group(fin,'MAIN') ) then
        call ffdev_utils_exit(DEV_OUT,1,'Specified control file does not contain {MAIN} group!')
    end if
    
    if( command_argument_count() .gt. 1 ) then
        call get_command_argument(2, OptTopName)
    end if
    
    if( command_argument_count() .gt. 2 ) then
        call get_command_argument(3, OptCrdName)
        OptRstName = OptCrdName
    end if 
    
    if( command_argument_count() .gt. 3 ) then
        call get_command_argument(4, OptRstName)
    end if       

    call ffdev_geoopt_ctrl_files(fin)
    call ffdev_geoopt_ctrl_minimize(fin)

    ! check if everything was read
    if( prmfile_count_ulines(fin) .ne. 0 ) then
        write(DEV_OUT,*)
        call prmfile_dump(fin,DEV_OUT,.true.)
        call ffdev_utils_exit(DEV_OUT,1,'Unprocessed lines found in the control file!')
    end if

    ! release the file
    call prmfile_clear(fin)

    ! load files and setup optimizer -------------------------------------------

    write(DEV_OUT,*)
    call ffdev_utils_heading(DEV_OUT,'Loading input data', ':')

    ! load topology
    write(DEV_OUT,*)
    call ffdev_utils_heading(DEV_OUT,'Simplified Topology','=')

    call ffdev_topology_init(top)
    call ffdev_topology_load(top,OptTopName)
    call ffdev_topology_info(top)

    ! load coordinates
    write(DEV_OUT,*)
    call ffdev_utils_heading(DEV_OUT,'XYZ Coordinates','=')

    call ffdev_geometry_init(geo)
    call ffdev_geometry_load_point(geo,OptCrdName)
    call ffdev_geometry_info_input(geo)
    call ffdev_gradient_allocate(geo)

    ! check coordinates and topology
    call ffdev_geometry_check_z(top,geo)

    ! trajectory
    call ffdev_geoopt_opentraj(top)

    ! run optimizer ------------------------------------------------------------
    write(DEV_OUT,*)
    call ffdev_utils_heading(DEV_OUT,'Geometry Optimization',':')
    call ffdev_topology_finalize_setup(top)
    call ffdev_geoopt_run(DEV_OUT,top,geo)

    ! finalize all -------------------------------------------------------------
    write(DEV_OUT,*)
    call ffdev_utils_heading(DEV_OUT,'Finalization',':')
    call ffdev_geometry_save_xyz(geo,OptRstName)

    ! close trajectory
    call ffdev_geoopt_closetraj

    call ffdev_utils_footer('MM Optimize')

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

10 format('    mmoptimize <ctrlfile> [topology [input [output]]]')

end subroutine print_usage

!===============================================================================

end program ffdev_mmoptimize


