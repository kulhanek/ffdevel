! ==============================================================================
! This file is part of FFDevel.
!    Copyright (C) 2018 Petr Kulhanek, kulhanek@chemi.muni.cz
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

module ffdev_mmd3_dat

use dftd3_api
use ffdev_constants

! access to dftd3 data
TYPE(dftd3_calc)    :: loc_dftd3_calc
TYPE(dftd3_input)   :: loc_dftd3_input

! damping modes
integer,parameter   :: MMD3_BM  = 1                 ! Becke–Johnson damping function (BJ)

! mmd3 setup
integer             :: mmd3_damping = MMD3_BM

! calculation of fractional coordination numbers
logical             :: mmd3_use_frac_cn = .true.    ! use fractional CN derived from bond_r0
real(DEVDP)         :: mmd3_k1 = 16.0               ! to derive fractional CN

! BJ damping parameters
real(DEVDP)         :: mmd3_bj_a1                   ! damping function parameters
real(DEVDP)         :: mmd3_bj_a2

! ==============================================================================

end module ffdev_mmd3_dat
