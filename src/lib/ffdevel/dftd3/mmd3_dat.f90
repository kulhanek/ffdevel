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

! calculation of fractional coordination numbers
logical             :: mmd3_use_frac_cn = .true.    ! use fractional CN derived from bond_r0
real(DEVDP)         :: mmd3_k1 = 16.0               ! to derive fractional CN

! MMD3 data ---------------------------------------------------------------------

! MMD3 pair, in internal units
type MMD3_PAIR_TYPE
    real(DEVDP)         :: c6ave
    real(DEVDP)         :: c6sig
    real(DEVDP)         :: c8ave
    real(DEVDP)         :: c8sig
    real(DEVDP)         :: Rc       ! critical radius for BJ
    integer             :: num
end type MMD3_PAIR_TYPE

! LJ pair, in internal units
type LJ_PAIR_TYPE
    real(DEVDP)         :: eps
    real(DEVDP)         :: r0
    real(DEVDP)         :: c6
    integer             :: num
end type LJ_PAIR_TYPE

logical                             :: mmd3_data_loaded = .false.
type(MMD3_PAIR_TYPE),allocatable    :: mmd3_pairs(:,:)   ! ntypes x ntypes
type(LJ_PAIR_TYPE),allocatable      :: fflj_pairs(:,:)   ! ntypes x ntypes

! ==============================================================================

end module ffdev_mmd3_dat
