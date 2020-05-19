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

program ffdev_compenetop_program

    use ffdev_sizes
    use ffdev_utils
    use ffdev_constants
    use ffdev_variables
    use ffdev_topology
    use ffdev_geometry
    use ffdev_energy

    implicit none
    character(len=MAX_PATH) :: topname1     ! input topology name #1
    character(len=MAX_PATH) :: topname2     ! input topology name #2
    character(len=MAX_PATH) :: crdname      ! input coordinate name
    type(TOPOLOGY)          :: top1
    type(TOPOLOGY)          :: top2
    type(GEOMETRY)          :: geo1
    type(GEOMETRY)          :: geo2
    ! --------------------------------------------------------------------------

    call ffdev_utils_header('Compare Energy between Two Topologies')

    ! test number of input arguments
    if( command_argument_count() .ne. 3 ) then
        call print_usage()
        call ffdev_utils_exit(DEV_OUT,1,'Incorrect number of arguments was specified (three expected)!')
    end if

    call get_command_argument(1, topname1)
    call get_command_argument(2, topname2)
    call get_command_argument(3, crdname)

    ! paths
    write(DEV_OUT,*)
    write(DEV_OUT,100) trim(topname1)
    write(DEV_OUT,110) trim(topname2)
    write(DEV_OUT,120) trim(crdname)

    ! load topology
    write(DEV_OUT,*)
    call ffdev_utils_heading(DEV_OUT,'Simplified Topology #1','=')

    call ffdev_topology_init(top1)
    call ffdev_topology_load(top1,topname1)
    call ffdev_topology_info(top1)

    ! load topology
    write(DEV_OUT,*)
    call ffdev_utils_heading(DEV_OUT,'Simplified Topology #2','=')

    call ffdev_topology_init(top2)
    call ffdev_topology_load(top2,topname2)
    call ffdev_topology_info(top2)

    ! load coordinates #1
    write(DEV_OUT,*)
    call ffdev_utils_heading(DEV_OUT,'Input coordinates','=')

    call ffdev_geometry_init(geo1)
    call ffdev_geometry_load_xyz(geo1,crdname)
    call ffdev_geometry_info_input(geo1)

    ! check coordinates and topology
    call ffdev_geometry_check_z(top1,geo1)
    call ffdev_geometry_check_z(top2,geo1)

    ! copy coordinates
    call ffdev_geometry_init(geo2)
    call ffdev_geometry_copy(geo2,geo1)

    ! finalize topology setup and calculate energy
    call ffdev_topology_finalize_setup(top1)
    call ffdev_topology_finalize_setup(top2)
    call ffdev_energy_all(top1,geo1)
    call ffdev_energy_all(top2,geo2)

    write(DEV_OUT,*)
    call ffdev_utils_heading(DEV_OUT,'FF Energies #1','=')
    call ffdev_geometry_info_ene(geo1)

    write(DEV_OUT,*)
    call ffdev_utils_heading(DEV_OUT,'FF Energies #2','=')
    call ffdev_geometry_info_ene(geo2)

    ! compare energy items
    write(DEV_OUT,*)
    call ffdev_utils_heading(DEV_OUT,'ENERGY','=')

    write(DEV_OUT,*)
    write(DEV_OUT,200)
    write(DEV_OUT,210)
    write(DEV_OUT,300) geo1%bond_ene,geo1%bond_ene,geo2%bond_ene-geo1%bond_ene
    write(DEV_OUT,310) geo1%angle_ene,geo2%angle_ene,geo2%angle_ene-geo1%angle_ene
    write(DEV_OUT,320) geo1%dih_ene,geo2%dih_ene,geo2%dih_ene-geo1%dih_ene
    write(DEV_OUT,330) geo1%impropr_ene,geo2%impropr_ene,geo2%impropr_ene-geo1%impropr_ene
    write(DEV_OUT,340) geo1%ele14_ene,geo2%ele14_ene,geo2%ele14_ene-geo1%ele14_ene
    write(DEV_OUT,350) geo1%nb14_ene,geo2%nb14_ene,geo2%nb14_ene-geo1%nb14_ene
    write(DEV_OUT,360) geo1%ele_ene,geo2%ele_ene,geo2%ele_ene-geo1%ele_ene
    write(DEV_OUT,370) geo1%nb_ene,geo2%nb_ene,geo2%nb_ene-geo1%nb_ene
    write(DEV_OUT,380) geo1%total_ene,geo2%total_ene,geo2%total_ene-geo1%total_ene


    call ffdev_utils_footer('Compare Energy between Two Topologies')

100 format('Simplified topology #1 : ',A)
110 format('Simplified topology #2 : ',A)
120 format('Input coordinates      : ',A)

200 format('#                        #1                   #2               diff(2-1)       ')
210 format('# ------------- -------------------- -------------------- -------------------- ')
300 format('# Ebonds     = ',F20.7,1X,F20.7,1X,F20.7)
310 format('# Eangles    = ',F20.7,1X,F20.7,1X,F20.7)
320 format('# Edihedrals = ',F20.7,1X,F20.7,1X,F20.7)
330 format('# Eimpropers = ',F20.7,1X,F20.7,1X,F20.7)
340 format('# E14ele     = ',F20.7,1X,F20.7,1X,F20.7)
350 format('# E14vdw     = ',F20.7,1X,F20.7,1X,F20.7)
360 format('# Eele       = ',F20.7,1X,F20.7,1X,F20.7)
370 format('# Evdw       = ',F20.7,1X,F20.7,1X,F20.7)
380 format('# Etotal     = ',F20.7,1X,F20.7,1X,F20.7)

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

10 format('    compenetop <stopology#1> <stopology#2> <xyzcrd>')

end subroutine print_usage

!===============================================================================

end program ffdev_compenetop_program
