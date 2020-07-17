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
! density overlap

integer,parameter   :: DO_PBE0_def2QZVPP = 1

integer             :: densoverlap_source            = DO_PBE0_def2QZVPP
logical             :: densoverlap_mod_by_charge     = .false.

! overlap of electron densities of isolated atoms ------------------------------

integer,parameter   :: DENSOVERLAP_MAX_Z = 86

! log(S_{rho}) = -b*r + a
! bp, ap - cation (+1)
! b0, a0 - neutral atom
! bm, am - anion (-1)

real(DEVDP)     :: densoverlap_bpii(1:DENSOVERLAP_MAX_Z)
real(DEVDP)     :: densoverlap_b0ii(1:DENSOVERLAP_MAX_Z)
real(DEVDP)     :: densoverlap_bmii(1:DENSOVERLAP_MAX_Z)

real(DEVDP)     :: densoverlap_apii(1:DENSOVERLAP_MAX_Z)
real(DEVDP)     :: densoverlap_a0ii(1:DENSOVERLAP_MAX_Z)
real(DEVDP)     :: densoverlap_amii(1:DENSOVERLAP_MAX_Z)

real(DEVDP)     :: densoverlap_b0ij(1:DENSOVERLAP_MAX_Z,1:DENSOVERLAP_MAX_Z)
real(DEVDP)     :: densoverlap_a0ij(1:DENSOVERLAP_MAX_Z,1:DENSOVERLAP_MAX_Z)



! ------------------------------------------------------------------------------

end module ffdev_densoverlap_dat
