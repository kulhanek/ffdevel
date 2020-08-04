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

module ffdev_atomicdata_dat

use ffdev_constants
use ffdev_variables

! ------------------------------------------------------------------------------
! density overlap - source data methods
integer,parameter   :: AO_PBE0_def2QZVPP    = 1

! modificators
integer,parameter   :: AO_MODS_PLAIN        = 1
integer,parameter   :: AO_MODS_BY_XDM       = 2

! ------------------------------------------------------------------------------
! setup

integer             :: atomoverlap_source   = AO_PBE0_def2QZVPP
integer             :: atomoverlap_mods     = AO_MODS_PLAIN

! ------------------------------------------------------------------------------
! databases

! overlap of electron densities of isolated atoms ------------------------------

! model log(Sij) =  a0 + log(1.0+b0*R+(b0*R)**2/3.0) - b0*R

integer,parameter   :: DENSOVERLAP_MAX_Z = 86

real(DEVDP)     :: densoverlap_b0ii(1:DENSOVERLAP_MAX_Z)
real(DEVDP)     :: densoverlap_a0ii(1:DENSOVERLAP_MAX_Z)

! overlap of WFs of isolated atoms ---------------------------------------------

! model log(Sij^2/R) = a0 + 2.0*log(1.0+b0*R/2.0+(b0*R)**2/12.0) -log(R)  - b0*R

integer,parameter   :: WFOVERLAP_MAX_Z = 86

real(DEVDP)     :: wfoverlap_b0ii(1:WFOVERLAP_MAX_Z)
real(DEVDP)     :: wfoverlap_a0ii(1:WFOVERLAP_MAX_Z)

! ------------------------------------------------------------------------------

end module ffdev_atomicdata_dat
