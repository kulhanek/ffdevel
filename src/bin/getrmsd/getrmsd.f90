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

program ffdev_getrmsd_program

    use ffdev_sizes
    use ffdev_utils
    use ffdev_constants
    use ffdev_geometry
    use ffdev_geometry_utils

    implicit none
    character(len=MAX_PATH) :: refname  ! reference coordinates
    character(len=MAX_PATH) :: srcname  ! reference coordinates
    type(GEOMETRY)          :: ref
    type(GEOMETRY)          :: src
    integer                 :: i
    ! --------------------------------------------------------------------------

    call ffdev_utils_header('Get RMSD')

    ! test number of input arguments
    if( (command_argument_count() .ne. 2) .and. (command_argument_count() .ne. 3)  ) then
        call print_usage()
        call ffdev_utils_exit(DEV_OUT,1,'Incorrect number of arguments was specified (two or three expected)!')
    end if

    call get_command_argument(1, refname)
    call get_command_argument(2, srcname)

    ! process input file -------------------------------------------------------
    write(DEV_OUT,*)
    call ffdev_utils_heading(DEV_OUT,'Reading coordinates files', ':')
    write(DEV_OUT,100) trim(refname)

    call ffdev_geometry_init(ref)
    call ffdev_geometry_load_point(ref,refname)

    write(DEV_OUT,110) trim(srcname)

    call ffdev_geometry_init(src)
    call ffdev_geometry_load_point(src,srcname)

    ! print input data ---------------------------------------------------------
    write(DEV_OUT,*)
    call ffdev_utils_heading(DEV_OUT,'Reference coordinates', ':')
    write(DEV_OUT,*)
    call ffdev_geometry_print_xyz(ref)

    write(DEV_OUT,*)
    call ffdev_utils_heading(DEV_OUT,'Source coordinates', ':')
    write(DEV_OUT,*)
    call ffdev_geometry_print_xyz(src)

    ! check compability of two geometries
    if( ref%natoms .ne. src%natoms ) then
        call ffdev_utils_exit(DEV_OUT,1,'Two geometries must have the same number of atoms!')
    end if

    do i=1,ref%natoms
        if( ref%z(i) .ne. src%z(i) ) then
            call ffdev_utils_exit(DEV_OUT,1,'Two geometries must have the same order of atoms!')
        end if
    end do

    ! results ------------------------------------------------------------------
    write(DEV_OUT,*)
    call ffdev_utils_heading(DEV_OUT,'Results', ':')
    write(DEV_OUT,*)
    if( command_argument_count() .eq. 2 ) then
        write(DEV_OUT,120) ffdev_geometry_utils_get_rmsd(ref%natoms,ref%z,ref%crd,src%crd,.true.)
    else
        write(DEV_OUT,130) ffdev_geometry_utils_get_rmsd_nofit(ref%natoms,ref%z,ref%crd,src%crd,.true.)
    end if

    ! end
    call ffdev_utils_footer('Get RMSD')

100 format('Reference coordinates : ',A)
110 format('Source coordinates    : ',A)
120 format('RMSD (fitted)                   = ',F16.6)
130 format('RMSD (nofit)                    = ',F16.6)

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

10 format('    getrmsd <reference> <source>')

end subroutine print_usage

!===============================================================================

end program ffdev_getrmsd_program
