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

program ffdev_gradient_program

    use ffdev_sizes
    use ffdev_utils
    use ffdev_constants
    use ffdev_topology
    use ffdev_geometry
    use ffdev_gradient
    use ffdev_gradient_utils

    implicit none

    character(len=MAX_PATH) :: topname      ! input topology name
    character(len=MAX_PATH) :: crdname      ! input coordinate name
    character(len=MAX_PATH) :: arg
    type(TOPOLOGY)          :: top
    type(GEOMETRY)          :: geo,ngeo
    logical                 :: do_test, do_numerical
    integer                 :: i
    ! --------------------------------------------------------------------------

    call ffdev_utils_header('RST Gradient')

    ! test number of input arguments
    if( command_argument_count() .lt. 2  ) then
        call print_usage()
        call ffdev_utils_exit(DEV_OUT,1,'Incorrect number of arguments was specified (at least two expected)!')
    end if

    call get_command_argument(1, topname)
    call get_command_argument(2, crdname)

    do_test = .false.
    do_numerical = .false.

    do i=3,command_argument_count()
        call get_command_argument(i, arg)
        select case(trim(arg))
            case('test')
                do_test = .true.
            case('numerical')
                do_numerical = .true.
            case default
                call ffdev_utils_exit(DEV_OUT,1,'Unrecognized argument ('//trim(arg)//')')
        end select
    end do

    ! paths
    write(DEV_OUT,*)
    write(DEV_OUT,100) trim(topname)
    write(DEV_OUT,110) trim(crdname)

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
    call ffdev_geometry_load_point(geo,crdname)
    call ffdev_geometry_info_input(geo)

    ! check coordinates and topology
    call ffdev_geometry_check_z(top,geo)

    ! finalize topology and geometry setup
    call ffdev_topology_finalize_setup(top)
    call ffdev_gradient_allocate(geo)
    call ffdev_geometry_init(ngeo)
    call ffdev_geometry_copy(ngeo,geo)
    call ffdev_gradient_allocate(ngeo)

    ! calculate energy and gradient
    ! only restraints if any
    write(DEV_OUT,*)
    if( do_numerical ) then
        write(DEV_OUT,'(A)') 'Numerical gradient ...'
        call ffdev_geometry_get_rst_penalty_num(geo)
    else
        write(DEV_OUT,'(A)') 'Analytical gradient ...'
        call ffdev_geometry_get_rst_penalty(geo)
    end if

    write(DEV_OUT,*)
    call ffdev_geometry_rstsum(DEV_OUT,geo)

    if( do_test ) then
        write(DEV_OUT,*)
        write(DEV_OUT,'(A)') '>>> INFO: Testing gradient ...'
        write(DEV_OUT,*)
        write(DEV_OUT,'(A)') 'Numerical gradient ...'
        call ffdev_geometry_get_rst_penalty_num(ngeo)

        if( .not. ffdev_gradient_test(geo,ngeo,1.0d-3) ) then
            write(DEV_OUT,*)
            call ffdev_utils_heading(DEV_OUT,'Analytical gradient','=')
            call ffdev_gradient_print(DEV_OUT,top,geo)
                write(DEV_OUT,*)
            call ffdev_utils_heading(DEV_OUT,'Numerical gradient','=')
            call ffdev_gradient_print(DEV_OUT,top,ngeo)

            write(DEV_OUT,*)
            call ffdev_utils_heading(DEV_OUT,'Difference in gradients','=')
            ngeo%grd = ngeo%grd - geo%grd
            call ffdev_gradient_print(DEV_OUT,top,ngeo)

            call ffdev_utils_exit(DEV_OUT,1,'Analytical and numerical gradients do not match!')
        else
            write(DEV_OUT,'(A)') 'Analytical and numerical gradients match each other ...'
            write(DEV_OUT,*)
            call ffdev_utils_heading(DEV_OUT,'Numerical gradient','=')
            call ffdev_gradient_print(DEV_OUT,top,ngeo)
        end if
    end if

    write(DEV_OUT,*)
    call ffdev_utils_heading(DEV_OUT,'Gradient','=')
    call ffdev_gradient_print(DEV_OUT,top,geo)

    call ffdev_utils_footer('RST Gradient')

100 format('Simplified topology : ',A)
110 format('Input point         : ',A)

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

10 format('    rstgradient <stopology> <pts> [test] [numerica]')

end subroutine print_usage

!===============================================================================

end program ffdev_gradient_program
