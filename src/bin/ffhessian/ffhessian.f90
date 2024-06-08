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

program ffdev_hessian_program

    use ffdev_sizes
    use ffdev_utils
    use ffdev_constants
    use ffdev_variables
    use ffdev_topology
    use ffdev_geometry
    use ffdev_energy
    use ffdev_gradient
    use ffdev_gradient_utils
    use ffdev_hessian
    use ffdev_hessian_utils
    use ffdev_timers

    implicit none

    character(len=MAX_PATH) :: topname      ! input topology name
    character(len=MAX_PATH) :: crdname      ! input coordinate name
    character(len=MAX_PATH) :: arg
    type(TOPOLOGY)          :: top
    type(GEOMETRY)          :: geo,ngeo,nggeo
    logical                 :: do_test,do_numerical,do_test_ene
    logical                 :: write_pts
    integer                 :: i
    real(DEVDP)             :: serr2
    ! --------------------------------------------------------------------------

    call ffdev_utils_header('FF Hessian')

    call ffdev_timers_init_top
    call ffdev_timers_init

    ! test number of input arguments
    if( command_argument_count() .lt. 2 ) then
        call print_usage()
        call ffdev_utils_exit(DEV_ERR,1,'Incorrect number of arguments was specified (at least two expected)!')
    end if

    call get_command_argument(1, topname)
    call get_command_argument(2, crdname)

    do_test = .false.
    do_test_ene = .false.
    do_numerical = .false.
    write_pts = .false.

    do i=3,command_argument_count()
        call get_command_argument(i, arg)
        select case(trim(arg))
            case('test')
                do_test = .true.
                do_test_ene = .false.
            case('test-viaene')
                do_test = .true.
                do_test_ene = .true.
            case('numerical')
                do_numerical = .true.
            case('write')
                write_pts = .true.
            case default
                call ffdev_utils_exit(DEV_ERR,1,'Unrecognized argument ('//trim(arg)//')')
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
    call ffdev_utils_heading(DEV_OUT,'XYZ Coordinates','=')

    call ffdev_geometry_init(geo)
    call ffdev_geometry_load_xyz(geo,crdname)
    call ffdev_geometry_info_input(geo)

    ! check coordinates and topology
    call ffdev_geometry_check_z(top,geo)

    ! finalize topology and geometry setup
    call ffdev_topology_finalize_setup(top)
    call ffdev_gradient_allocate(geo)
    call ffdev_hessian_allocate(geo)

    ! calculate energy and gradient
    write(DEV_OUT,*)
    if( do_numerical ) then
        write(DEV_OUT,'(A)') 'Numerical hessian ...'
        call ffdev_hessian_num_by_grds_all(top,geo)
    else
        write(DEV_OUT,'(A)') 'Analytical hessian ...'
        call ffdev_hessian_all(top,geo)
    end if

    if( write_pts ) then
        write(DEV_OUT,*)
        call ffdev_utils_heading(DEV_OUT,'Writing final data','=')
        geo%title = 'calculated by ffhessian from '//trim(topname)//' and '//trim(crdname)
        call ffdev_geometry_info_input(geo)
        call ffdev_geometry_save_point(geo,crdname,.false.)
    end if

    if( do_test ) then
        write(DEV_OUT,*)
        write(DEV_OUT,'(A)') '>>> INFO: Testing Hessian and gradient ...'
        call ffdev_geometry_init(nggeo)
        call ffdev_geometry_copy(nggeo,geo)
        call ffdev_gradient_allocate(nggeo)

        ! calculate analytical gradient
        write(DEV_OUT,*)
        write(DEV_OUT,'(A)') 'Analytical gradient from gradient calculation ...'
        call ffdev_gradient_all(top,nggeo)

        if( .not. ffdev_gradient_test(geo,nggeo,1.0d-3) ) then
            write(DEV_OUT,*)
            call ffdev_utils_heading(DEV_OUT,'Analytical FF Gradient from Hessian calculation','=')
            call ffdev_gradient_print(DEV_OUT,top,geo)
                write(DEV_OUT,*)
            call ffdev_utils_heading(DEV_OUT,'Analytical FF Gradient from gradient calculation','=')
            call ffdev_gradient_print(DEV_OUT,top,nggeo)

            write(DEV_OUT,*)
            call ffdev_utils_heading(DEV_OUT,'Difference in Gradients','=')
            nggeo%grd = nggeo%grd - geo%grd
            call ffdev_gradient_print(DEV_OUT,top,nggeo)

            call ffdev_utils_exit(DEV_ERR,1,'Analytical gradients from Hessian and gradient calculations do not match!')
        else
            write(DEV_OUT,'(A)') 'Analytical gradients from Hessian and gradient calculations match each other ...'
            write(DEV_OUT,*)
            call ffdev_utils_heading(DEV_OUT,'Analytical FF Gradient from gradient calculation','=')
            call ffdev_gradient_print(DEV_OUT,top,nggeo)
            write(DEV_OUT,*)
        end if

        call ffdev_geometry_init(ngeo)
        call ffdev_geometry_copy(ngeo,geo)
        call ffdev_gradient_allocate(ngeo)
        call ffdev_hessian_allocate(ngeo)

        write(DEV_OUT,'(A)') 'Numerical hessian ...'
        if( do_test_ene ) then
            write(DEV_OUT,'(A)') '> Employing energies only ...'
            call ffdev_hessian_num_all(top,ngeo)
        else
            write(DEV_OUT,'(A)') '> Employing gradients ...'
            call ffdev_hessian_num_by_grds_all(top,ngeo)
        end if

        if( .not. ffdev_hessian_test(geo,ngeo,1.0d-3) ) then
            write(DEV_OUT,*)
            call ffdev_utils_heading(DEV_OUT,'Analytical FF Hessian','=')
            call ffdev_hessian_print(DEV_OUT,geo)
            write(DEV_OUT,*)
            call ffdev_utils_heading(DEV_OUT,'Numerical FF Hessian','=')
            call ffdev_hessian_print(DEV_OUT,ngeo)
            write(DEV_OUT,*)
            call ffdev_utils_heading(DEV_OUT,'Difference in Hessians','=')
            ngeo%hess = ngeo%hess - geo%hess
            call ffdev_hessian_print(DEV_OUT,ngeo)

            serr2 = sum( ngeo%hess**2 )
            write(DEV_ERR,*) 'SERR2 = ', serr2

            call ffdev_utils_exit(DEV_ERR,1,'Analytical and numerical Hessians do not match!')
        else
            write(DEV_OUT,'(A)') 'Analytical and numerical Hessians match each other ...'
            write(DEV_OUT,*)
            call ffdev_utils_heading(DEV_OUT,'Numerical FF Hessian','=')
            call ffdev_hessian_print(DEV_OUT,ngeo)
        end if
    end if

    write(DEV_OUT,*)
    call ffdev_utils_heading(DEV_OUT,'FF Energies','=')
    call ffdev_geometry_info_ene(geo)

    write(DEV_OUT,*)
    call ffdev_utils_heading(DEV_OUT,'FF Gradient','=')
    call ffdev_gradient_print(DEV_OUT,top,geo)

    write(DEV_OUT,*)
    call ffdev_utils_heading(DEV_OUT,'FF Hessian','=')
    call ffdev_hessian_print(DEV_OUT,geo)

    write(DEV_OUT,*)
    call ffdev_utils_heading(DEV_OUT,'Vibrations','=')
    call ffdev_hessian_allocate_freq(geo)
    call ffdev_hessian_calc_freqs(geo)
    call ffdev_hessian_print_freqs(DEV_OUT,geo)
    write(DEV_OUT,*)
    call ffdev_hessian_print_freqs_lin(DEV_OUT,geo)

    if( ffdev_gradient_rmsg_only(geo) .gt. 0.1d0 ) then
        write(DEV_OUT,*)
        write(DEV_OUT,120)
    end if

    call ffdev_timers_finalize(.true.)

    call ffdev_utils_footer('FF Hessian')

100 format('Simplified topology : ',A)
110 format('Input coordinates   : ',A)
120 format('  !!!!! WARNING !!!!! Frequencies might be inaccurate due to high gradient.')

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

10 format('    ffhessian <stopology> <point> [test] [write] [numerical]')

end subroutine print_usage

!===============================================================================

end program ffdev_hessian_program
