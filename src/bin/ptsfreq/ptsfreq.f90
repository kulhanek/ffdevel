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

program ffdev_ptsfreq_program

    use ffdev_sizes
    use ffdev_utils
    use ffdev_constants
    use ffdev_variables
    use ffdev_geometry
    use ffdev_hessian_utils
    use ffdev_gradient_utils

    implicit none
    character(len=MAX_PATH) :: geoname  ! input point file name
    type(GEOMETRY)          :: geo
    ! --------------------------------------------------------------------------

    call ffdev_utils_header('Point Frequencies')

    ! test number of input arguments
    if( command_argument_count() .ne. 1 ) then
        call print_usage()
        call ffdev_utils_exit(DEV_ERR,1,'Incorrect number of arguments was specified (one expected)!')
    end if

    call get_command_argument(1, geoname)

    ! process input file -------------------------------------------------------
    write(DEV_OUT,*)
    call ffdev_utils_heading(DEV_OUT,'Reading point file', ':')
    write(DEV_OUT,'(a,a)') 'Point file name : ',trim(geoname)

    call ffdev_geometry_init(geo)
    call ffdev_geometry_load_point(geo,geoname)

    write(DEV_OUT,*)
    call ffdev_geometry_info_point_header(GEO_INFO_NOENERGY)
    call ffdev_geometry_info_point(geo,GEO_INFO_NOENERGY)

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
            write(DEV_OUT,120)
        end if
    end if

    ! end
    call ffdev_utils_footer('Point Frequencies')

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

10 format('    ptfreq <point>')

end subroutine print_usage

!===============================================================================

end program ffdev_ptsfreq_program
