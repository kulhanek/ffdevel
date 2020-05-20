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

program ffdev_compgeo_program

    use ffdev_sizes
    use ffdev_utils
    use ffdev_constants
    use ffdev_variables
    use ffdev_topology
    use ffdev_geometry
    use ffdev_geometry_utils
    use ffdev_energy

    implicit none
    character(len=MAX_PATH) :: topname      ! input topology name
    character(len=MAX_PATH) :: crdname1     ! input coordinate name #1
    character(len=MAX_PATH) :: crdname2     ! input coordinate name #2
    type(TOPOLOGY)          :: top
    type(GEOMETRY)          :: geo1
    type(GEOMETRY)          :: geo2
    ! --------------------------------------------------------------------------

    call ffdev_utils_header('Compare Geometry')

    ! test number of input arguments
    if( command_argument_count() .ne. 3 ) then
        call print_usage()
        call ffdev_utils_exit(DEV_ERR,1,'Incorrect number of arguments was specified (three expected)!')
    end if

    call get_command_argument(1, topname)
    call get_command_argument(2, crdname1)
    call get_command_argument(3, crdname2)

    ! paths
    write(DEV_OUT,*)
    write(DEV_OUT,100) trim(topname)
    write(DEV_OUT,110) trim(crdname1)
    write(DEV_OUT,120) trim(crdname2)

    ! load topology
    write(DEV_OUT,*)
    call ffdev_utils_heading(DEV_OUT,'Simplified Topology','=')

    call ffdev_topology_init(top)
    call ffdev_topology_load(top,topname)
    call ffdev_topology_info(top)

    ! load coordinates #1
    write(DEV_OUT,*)
    call ffdev_utils_heading(DEV_OUT,'XYZ Coordinates #1','=')

    call ffdev_geometry_init(geo1)
    call ffdev_geometry_load_xyz(geo1,crdname1)
    call ffdev_geometry_info_input(geo1)

    ! check coordinates and topology
    call ffdev_geometry_check_z(top,geo1)

    ! load coordinates #2
    write(DEV_OUT,*)
    call ffdev_utils_heading(DEV_OUT,'XYZ Coordinates #2','=')

    call ffdev_geometry_init(geo2)
    call ffdev_geometry_load_xyz(geo2,crdname2)
    call ffdev_geometry_info_input(geo2)

    ! check coordinates and topology
    call ffdev_geometry_check_z(top,geo2)

    ! finalize topology setup and calculate energy
    call ffdev_topology_finalize_setup(top)

    ! compare geometry
    write(DEV_OUT,*)
    call ffdev_utils_heading(DEV_OUT,'BONDS','=')
    call ffdev_geometry_utils_comp_bonds(.true.,top,geo1%crd,geo2%crd)
    write(DEV_OUT,*)
    call ffdev_utils_heading(DEV_OUT,'ANGLES','=')
    call ffdev_geometry_utils_comp_angles(.true.,top,geo1%crd,geo2%crd)
    write(DEV_OUT,*)
    call ffdev_utils_heading(DEV_OUT,'DIHEDRALS','=')
    call ffdev_geometry_utils_comp_dihedrals(.true.,top,geo1%crd,geo2%crd,.false.)

    call ffdev_utils_footer('Compare Geometry')

100 format('Simplified topology : ',A)
110 format('XYZ coordinates #1  : ',A)
120 format('XYZ coordinates #2  : ',A)

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

10 format('    compgeo <stopology> <xyzcrd#1> <xyzcrd#2>')

end subroutine print_usage

!===============================================================================

end program ffdev_compgeo_program
