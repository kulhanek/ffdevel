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

integer,parameter   :: PAULI_MAX_NSTO  = 3          ! max allowed number of STO orbitals

! types of WF
integer,parameter   :: PAULI_WF_SV     = 1          ! split valence 1s STO orbitals
integer,parameter   :: PAULI_WF_STO    = 2          ! STO type orbitals 1s,2s,3s,...

! types of XFUN
integer,parameter   :: PAULI_XFUN_RHOP = 1          ! power of density
integer,parameter   :: PAULI_XFUN_KXEG  = 2         ! kinetic-exchange for electron gas

! numgrid setup
logical     :: pauli_use_numgrid    ! use numgrid for integration
logical     :: pauli_cache_grid     ! cache grid for integration

! WFN
integer     :: pauli_wf_nsto        ! number of STO orbitals
logical     :: pauli_wf_truncbyn    ! truncate nsto by n (main quantum number)
integer     :: pauli_wf2rho_power   ! wf->rho transformation
integer     :: pauli_wf_form        ! type of WF

! PDENS
real(DEVDP) :: pauli_dens_dpower    ! can be optimized via pauli_dp
real(DEVDP) :: pauli_dens_rpower    ! can be optimized via pauli_rp

! XFUN
integer     :: pauli_xfun           ! type of xr functional
real(DEVDP) :: pauli_xfun_dpower    ! can be optimized via pauli_xd
real(DEVDP) :: pauli_xfun_kpower    ! can be optimized via pauli_xk
real(DEVDP) :: pauli_xfun_xpower    ! can be optimized via pauli_xx
real(DEVDP) :: pauli_xfun_xfac      ! can be optimized via pauli_xf

end module ffdev_pauli_dat
