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

module ffdev_dftd3

use dftd3_api
use ffdev_constants

logical             :: dftd3_initialized = .false.
TYPE(dftd3_calc)    :: loc_dftd3_calc
TYPE(dftd3_input)   :: loc_dftd3_input

contains

! ==============================================================================
! subroutine ffdev_dftd3_get_c6
! ==============================================================================

real(DEVDP) function ffdev_dftd3_get_c6(z1,cn1,z2,cn2)

    implicit none
    integer         :: z1,z2
    real(DEVDP)     :: cn1,cn2
    ! --------------------------------------------
    real(DEVDP)     :: c6
    ! --------------------------------------------------------------------------

    if( .not. dftd3_initialized ) then
        call dftd3_init(loc_dftd3_calc,loc_dftd3_input)
        dftd3_initialized = .true.
    end if

    call getc6(maxc,max_elem,loc_dftd3_calc%c6ab,loc_dftd3_calc%mxc,z1,z2,cn1,cn2,c6)

    ffdev_dftd3_get_c6 = c6 * 627.50960803059 / 1.889725989d0**6

end function ffdev_dftd3_get_c6

! ==============================================================================
! subroutine ffdev_dftd3_get_c8
! ==============================================================================

real(DEVDP) function ffdev_dftd3_get_c8(z1,cn1,z2,cn2)

    implicit none
    integer         :: z1,z2
    real(DEVDP)     :: cn1,cn2
    ! --------------------------------------------
    real(DEVDP)     :: c6
    ! --------------------------------------------------------------------------

    if( .not. dftd3_initialized ) then
        call dftd3_init(loc_dftd3_calc,loc_dftd3_input)
        dftd3_initialized = .true.
    end if

    call getc6(maxc,max_elem,loc_dftd3_calc%c6ab,loc_dftd3_calc%mxc,z1,z2,cn1,cn2,c6)

    ffdev_dftd3_get_c8 = 3.0d0*c6*r2r4(z1)*r2r4(z2) * 627.50960803059 / 1.889725989d0**8

end function ffdev_dftd3_get_c8

! ==============================================================================

end module ffdev_dftd3
