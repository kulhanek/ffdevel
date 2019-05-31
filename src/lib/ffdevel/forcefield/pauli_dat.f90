! ==============================================================================
! This file is part of FFDevel.
!    Copyright (C) 2019 Petr Kulhanek, kulhanek@chemi.muni.cz
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

module ffdev_pauli_dat

use ffdev_sizes

! spli-valence modes
integer,parameter   :: PAULI_SZ   = 1
integer,parameter   :: PAULI_DZ   = 2

! Pauili repulsion
logical     :: pauli_use_numgrid    = .true.    ! use numgrid for integration
logical     :: pauli_cache_grid     = .true.    ! cache grid for Eexchrep calculation
real(DEVDP) :: pauli_dens_power     = 1.0       ! can be optimized via pauli_dp
real(DEVDP) :: pauli_lda_power      = 2.0       ! can be optimized via pauli_lp
integer     :: pauli_split_valence  = PAULI_DZ

end module ffdev_pauli_dat
