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

module ffdev_densoverlap_dat

use ffdev_constants
use ffdev_variables

! ------------------------------------------------------------------------------

integer,parameter   :: DO_HF_DKH_ANORCC = 1

integer             :: densoverlap_source        = DO_HF_DKH_ANORCC
logical             :: densoverlap_mod_by_charge = .false.

! atom densities ---------------------------------------------------------------

integer,parameter   :: DENSOVERLAP_MAX_Z = 86

! log(S_{rho}) = -b*r + a
! bp, ap - cation (+1)
! b0, a0 - neutral atom
! bm, am - anion (-1)

real(DEVDP)     :: densoverlap_bp(1:DENSOVERLAP_MAX_Z)
real(DEVDP)     :: densoverlap_b0(1:DENSOVERLAP_MAX_Z)
real(DEVDP)     :: densoverlap_bm(1:DENSOVERLAP_MAX_Z)

real(DEVDP)     :: densoverlap_ap(1:DENSOVERLAP_MAX_Z)
real(DEVDP)     :: densoverlap_a0(1:DENSOVERLAP_MAX_Z)
real(DEVDP)     :: densoverlap_am(1:DENSOVERLAP_MAX_Z)

! ------------------------------------------------------------------------------

end module ffdev_densoverlap_dat
