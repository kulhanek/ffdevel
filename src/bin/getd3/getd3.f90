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

program ffdev_getd3_program

    use ffdev_sizes
    use ffdev_utils
    use ffdev_constants
    use ffdev_mmd3
    use smf_periodic_table

    implicit none
    logical                     :: rst
    character(SMF_MAX_SYMBOL)   :: ele1
    character(SMF_MAX_SYMBOL)   :: ele2
    character(MAX_PATH)         :: tmp
    integer                     :: z1, z2
    real(DEVDP)                 :: c6, c8, c10, r86, r108, cn1, cn2
    ! --------------------------------------------------------------------------

    ! test number of input arguments
    if( command_argument_count() .ne. 4 ) then
        call print_usage()
        call ffdev_utils_exit(DEV_OUT,1,'Incorrect number of arguments was specified (four expected)!')
    end if

    call get_command_argument(1, ele1)
    write(DEV_OUT,100) trim(ele1)

    z1 = SearchZBySymbol(ele1)
    if( z1 .eq. 0 ) then
        call ffdev_utils_exit(DEV_OUT,1,'Unable to decode "'// trim(ele1) // '" element!')
    end if
    write(DEV_OUT,105) z1

    call get_command_argument(2, tmp)
    read(tmp,*) cn1
    write(DEV_OUT,106) cn1

    call get_command_argument(3, ele2)
    write(DEV_OUT,110) trim(ele2)

    z2 = SearchZBySymbol(ele2)
    if( z2 .eq. 0 ) then
        call ffdev_utils_exit(DEV_OUT,1,'Unable to decode "'// trim(ele2) // '" element!')
    end if
    write(DEV_OUT,115) z2

    call get_command_argument(4, tmp)
    read(tmp,*) cn2
    write(DEV_OUT,116) cn2

    ! init mmd3
    call ffdev_mmd3_init

    c6 = ffdev_mmd3_get_c6_by_z(z1,cn1,z2,cn2)
    c8 = ffdev_mmd3_get_c8_by_z(z1,cn1,z2,cn2)
    c10 = (49.0d0/40.0d0)*c8**2/c6
    r86 = sqrt(c8/c6)
    r108 = sqrt(c10/c8)

    write(DEV_OUT,120) c6
    write(DEV_OUT,130) c8
    write(DEV_OUT,140) c10
    write(DEV_OUT,150) r86
    write(DEV_OUT,160) r108

100 format('Element#1: ',A3)
105 format('Z#1:       ',I3)
106 format('CN#1:      ',F5.3)
110 format('Element#2: ',A3)
115 format('Z#2:       ',I3)
116 format('CN#2:      ',F5.3)
120 format('C6:        ',E24.16)
130 format('C8:        ',E24.16)
140 format('C10:       ',E24.16)
150 format('R8/6:      ',E24.16)
160 format('R10/8:     ',E24.16)

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

10 format('    getd3 <ele1> <cn1> <ele2> <cn2>')

end subroutine print_usage

!===============================================================================

end program ffdev_getd3_program
