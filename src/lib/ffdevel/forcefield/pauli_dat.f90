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

! types of WF
integer,parameter   :: PAULI_WF_SV      = 1     ! split valence an*exp(-bn*r)
integer,parameter   :: PAULI_WF_RSV     = 2     ! split valence an*exp(-bn*r)*r**cn

! types of XFUN
integer,parameter   :: PAULI_XFUN_RHOP  = 1     ! power of density
integer,parameter   :: PAULI_XFUN_KXEG  = 2     ! kinetic-exchange for electron gas

! numgrid setup
logical     :: pauli_use_numgrid    ! use numgrid for integration
logical     :: pauli_cache_grid     ! cache grid for integration

! WFN
integer     :: pauli_wf_form        ! type of WF
integer     :: pauli_wf_nsto        ! number of orbitals
integer     :: pauli_wf2rho_power   ! wf->rho transformation (from setup)

! PDENS
real(DEVDP) :: pauli_dens_dpower    ! can be optimized via pauli_dp
real(DEVDP) :: pauli_dens_rpower    ! can be optimized via pauli_rp

! XFUN
integer     :: pauli_xfun           ! type of xr functional
real(DEVDP) :: pauli_xfun_dpower    ! can be optimized via pauli_xd
real(DEVDP) :: pauli_xfun_kpower    ! can be optimized via pauli_xk
real(DEVDP) :: pauli_xfun_xpower    ! can be optimized via pauli_xx
real(DEVDP) :: pauli_xfun_xfac      ! can be optimized via pauli_xf

! working data
integer     :: used_wf2rho_power    ! actually used value, setup in preform_parameters

end module ffdev_pauli_dat
