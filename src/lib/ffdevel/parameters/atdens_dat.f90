! ==============================================================================
! This file is part of FFDevel.
!    Copyright (C) 2020 Petr Kulhanek, kulhanek@chemi.muni.cz
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

module ffdev_atdens_dat

use ffdev_constants
use ffdev_variables

! ------------------------------------------------------------------------------

integer,parameter   :: ATDENS_HF_UGBS       = 1
integer,parameter   :: ATDENS_CC_UGBS       = 2
integer,parameter   :: ATDENS_HF_DKH_ANORCC = 3
integer,parameter   :: ATDENS_CC_DKH_ANORCC = 4

integer             :: atdens_source        = ATDENS_HF_DKH_ANORCC
logical             :: atdens_mod_by_charge = .false.

! atom densities ---------------------------------------------------------------

integer,parameter   :: ATDENS_MAX_Z = 86

! log(rho) = -b*r + a
! bp, ap - cationt (+1)
! b0, a0 - neutral atom
! bm, am - anion (-1)

real(DEVDP)     :: atdens_bp(1:ATDENS_MAX_Z)
real(DEVDP)     :: atdens_b0(1:ATDENS_MAX_Z)
real(DEVDP)     :: atdens_bm(1:ATDENS_MAX_Z)

real(DEVDP)     :: atdens_ap(1:ATDENS_MAX_Z)
real(DEVDP)     :: atdens_a0(1:ATDENS_MAX_Z)
real(DEVDP)     :: atdens_am(1:ATDENS_MAX_Z)

! ------------------------------------------------------------------------------

end module ffdev_atdens_dat
