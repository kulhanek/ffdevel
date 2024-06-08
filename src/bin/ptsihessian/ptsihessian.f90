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

program ffdev_ptsihessian_program

    use ffdev_sizes
    use ffdev_utils
    use ffdev_constants
    use ffdev_variables
    use ffdev_topology
    use ffdev_geometry
    use ffdev_hessian_utils
    use ffdev_gradient_utils
    use ffdev_hessian
    use ffdev_nbmode_LJ

    implicit none
    character(len=MAX_PATH) :: topname
    character(len=MAX_PATH) :: geoname
    character(len=MAX_PATH) :: arg
    type(TOPOLOGY)          :: top
    type(GEOMETRY)          :: geo,geo_backup
    logical                 :: exclude_bonds,exclude_angles,exclude_dihedrals,exclude_impropers,exclude_nb
    integer                 :: i
    ! --------------------------------------------------------------------------

    call ffdev_utils_header('Hessian in Internal Coordinates')

    ! test number of input arguments
    if( command_argument_count() .lt. 2 ) then
        call print_usage()
        call ffdev_utils_exit(DEV_ERR,1,'Incorrect number of arguments was specified (at least two expected)!')
    end if

    call get_command_argument(1, topname)
    call get_command_argument(2, geoname)

    exclude_bonds = .false.
    exclude_angles = .false.
    exclude_dihedrals = .true.
    exclude_impropers = .true.
    exclude_nb = .true.

    do i=3,command_argument_count()
        call get_command_argument(i, arg)
        select case(trim(arg))
            case('exclude_bonds')
                exclude_bonds = .true.
            case('noexclude_angles')
                exclude_angles = .true.
            case('noexclude_dihedrals')
                exclude_dihedrals = .false.
            case('noexclude_impropers')
                exclude_impropers = .false.
            case('noexclude_nb')
                exclude_nb = .false.
            case default
                call ffdev_utils_exit(DEV_ERR,1,'Unrecognized argument ('//trim(arg)//')')
        end select
    end do

    ! paths
    write(DEV_OUT,*)
    write(DEV_OUT,100) trim(topname)
    write(DEV_OUT,110) trim(geoname)

    ! load topology
    write(DEV_OUT,*)
    call ffdev_utils_heading(DEV_OUT,'Simplified Topology','=')

    call ffdev_topology_init(top)
    call ffdev_topology_load(top,topname)
    call ffdev_topology_info(top)

    ! load coordinates
    write(DEV_OUT,*)
    call ffdev_utils_heading(DEV_OUT,'Input Coordinates','=')

    call ffdev_geometry_init(geo)
    call ffdev_geometry_load_point(geo,geoname)
    call ffdev_geometry_info_input(geo)

    ! check coordinates and topology
    call ffdev_geometry_check_z(top,geo)

    ! finalize topology setup and calculate energy
    call ffdev_topology_finalize_setup(top)

    if( .not. geo%trg_hess_loaded ) then
        call ffdev_utils_exit(DEV_ERR,1,'Point does not contain Hessian!')
    end if

    ! print input data ---------------------------------------------------------
    write(DEV_OUT,*)
    call ffdev_utils_heading(DEV_OUT,'Input coordinates', ':')

    write(DEV_OUT,*)
    call ffdev_geometry_print_xyz(geo)

    if( geo%trg_grd_loaded ) then
        ! allocate gradient and copy target one to it
        call ffdev_gradient_allocate(geo)
        geo%grd = geo%trg_grd

        write(DEV_OUT,*)
        call ffdev_utils_heading(DEV_OUT,'Input Gradient', ':')
        write(DEV_OUT,*)
        call ffdev_gradient_print_notop(DEV_OUT,geo)
    end if

    ! allocate hessian and copy target one to it
    call ffdev_hessian_allocate(geo)
    call ffdev_hessian_allocate_freq(geo)
    geo%hess = geo%trg_hess

    write(DEV_OUT,*)
    call ffdev_utils_heading(DEV_OUT,'Input Hessian', ':')
    write(DEV_OUT,*)
    call ffdev_hessian_print(DEV_OUT,geo)

    ! calculate frequencies ----------------------------------------------------
    call ffdev_hessian_calc_freqs(geo)

    ! print frequencies --------------------------------------------------------
    write(DEV_OUT,*)
    call ffdev_utils_heading(DEV_OUT,'Final Frequencies', ':')
    write(DEV_OUT,*)
    call ffdev_hessian_print_freqs(DEV_OUT,geo)
    write(DEV_OUT,*)
    call ffdev_hessian_print_freqs_lin(DEV_OUT,geo)

    if( geo%trg_grd_loaded ) then
        if( ffdev_gradient_rmsg_only(geo) .gt. 0.1d0 ) then
            write(DEV_OUT,*)
            write(DEV_OUT,130)
        end if
    end if

    ! hessian in internal coordinates
    write(DEV_OUT,*)
    call ffdev_utils_heading(DEV_OUT,'Internal Coordinates - Force Constants', ':')

    call ffdev_hessian_allocate_trg_ihess(top,geo)

    if( exclude_bonds .or. exclude_angles .or. exclude_dihedrals .or. exclude_impropers .or. exclude_nb ) then
        call ffdev_geometry_copy(geo_backup,geo)
        call ffdev_gradient_allocate(geo_backup)
        geo_backup%grd = 0.0d0
        call ffdev_hessian_allocate(geo_backup)
        geo_backup%hess = 0.0d0
    end if

    if( exclude_bonds ) then
        call ffdev_hessian_bonds(top,geo_backup)
        write(DEV_OUT,'(a)') '> Bonds excluded ...'
    end if

    if( exclude_angles ) then
        call ffdev_hessian_angles(top,geo_backup)
        write(DEV_OUT,'(a)') '> Angles excluded ...'
    end if

    if( exclude_dihedrals ) then
        call ffdev_hessian_dihedrals(top,geo_backup)
        write(DEV_OUT,'(a)') '> Dihedrals excluded ...'
    end if

    if( exclude_impropers ) then
        call ffdev_hessian_impropers(top,geo_backup)
        write(DEV_OUT,'(a)') '> Impropers excluded ...'
    end if

    if( exclude_nb ) then
        call ffdev_hessian_nb_lj(top,geo_backup)
        write(DEV_OUT,'(a)') '> NB excluded ...'
    end if

    if( exclude_bonds .or. exclude_angles .or. exclude_dihedrals .or. exclude_impropers .or. exclude_nb ) then
        geo%trg_ihess = geo%trg_hess - geo_backup%hess
        call ffdev_geometry_copy(geo,geo_backup)
    else
        geo%trg_ihess = geo%trg_hess
    end if

    ! calculate ihess
    call ffdev_hessian_calc_trg_ihess(top,geo)

    ! print results
    call ffdev_hessian_print_trg_ihess_bonds(top,geo)
    call ffdev_hessian_print_trg_ihess_angles(top,geo)

    ! end
    call ffdev_utils_footer('Hessian in Internal Coordinates')

100 format('Simplified topology : ',A)
110 format('Input coordinates   : ',A)

130 format('  !!!!! WARNING !!!!! Frequencies might be inaccurate due to high gradient.')

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

10 format('    ptsihessian <stop> <point> [excludenb]')

end subroutine print_usage

!===============================================================================

end program ffdev_ptsihessian_program
