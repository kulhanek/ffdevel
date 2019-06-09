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

module ffdev_xdm_dat

use ffdev_constants

! XDM data ---------------------------------------------------------------------

! XDM pair
type XDM_PAIR_TYPE
    real(DEVDP)         :: c6ave
    real(DEVDP)         :: c6sig
    real(DEVDP)         :: c8ave
    real(DEVDP)         :: c8sig
    real(DEVDP)         :: c10ave
    real(DEVDP)         :: c10sig
    integer             :: num
    real(DEVDP)         :: eps      ! LJ eps = C6/Rvdw**6
    real(DEVDP)         :: Rvdw     ! vdW radius derived from V, V0, and pol0
    real(DEVDP)         :: Cn(20)
end type XDM_PAIR_TYPE

! XDM atom
type XDM_ATOM_TYPE
    real(DEVDP)         :: vave     ! atom volume
    real(DEVDP)         :: vsig
    real(DEVDP)         :: v0ave    ! free atom volume
    real(DEVDP)         :: v0sig
    real(DEVDP)         :: p0ave    ! free atom polarizability
    real(DEVDP)         :: p0sig
    real(DEVDP)         :: pol      ! atomic polarizability
    integer             :: num
    real(DEVDP)         :: Rvdw     ! vdW radius derived from V, V0, and pol0
end type XDM_ATOM_TYPE

logical                         :: xdm_data_loaded = .false.
type(XDM_ATOM_TYPE),allocatable :: xdm_atoms(:)     ! ntypes
type(XDM_PAIR_TYPE),allocatable :: xdm_pairs(:,:)   ! ntypes x ntypes

! possible apply modes for LJ parameters from XDM data
integer,parameter       :: APPLY_XDM_NULL = 0
integer,parameter       :: APPLY_XDM_EPS  = 1
integer,parameter       :: APPLY_XDM_R0   = 2

! keep C6 constant modes
integer,parameter       :: LEFT_XDM_C6              = 0
integer,parameter       :: KEEP_XDM_C6_VIA_R0       = 1
integer,parameter       :: KEEP_XDM_C6_VIA_EPS      = 2

! === [xdm] ====================================================================
real(DEVDP)     :: xdm_rvdw_fac     =  2.54d0   ! FIXME -
real(DEVDP)     :: xdm_C6Scale      =  2.119    ! FIXME - LJ(C6)vsPotMinima
integer         :: xdm_C6Mode       = LEFT_XDM_C6

! ------------------------------------------------------------------------------

end module ffdev_xdm_dat
