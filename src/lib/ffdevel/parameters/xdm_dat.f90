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
use ffdev_variables

! XDM data ---------------------------------------------------------------------

! XDM pair, in internal units
type XDM_PAIR_TYPE
    real(DEVDP)         :: c6ave
    real(DEVDP)         :: c6sig
    real(DEVDP)         :: c8ave
    real(DEVDP)         :: c8sig
    real(DEVDP)         :: c10ave
    real(DEVDP)         :: c10sig
    real(DEVDP)         :: Rcd      ! critical from disp
    real(DEVDP)         :: Rcp      ! critical from pol (Rvdw)
    real(DEVDP)         :: Rc       ! critical radius for BJ
    integer             :: num
end type XDM_PAIR_TYPE

! XDM atom, in atomic units
type XDM_ATOM_TYPE
    real(DEVDP)         :: vave     ! atom volume
    real(DEVDP)         :: vsig
    real(DEVDP)         :: v0ave    ! free atom volume
    real(DEVDP)         :: v0sig
    real(DEVDP)         :: p0ave    ! free atom polarizability
    real(DEVDP)         :: p0sig
    real(DEVDP)         :: pol      ! atomic polarizability
    integer             :: num
end type XDM_ATOM_TYPE

logical                         :: xdm_data_loaded = .false.
type(XDM_ATOM_TYPE),allocatable :: xdm_atoms(:)     ! ntypes
type(XDM_PAIR_TYPE),allocatable :: xdm_pairs(:,:)   ! ntypes x ntypes

! ------------------------------------------------------------------------------

! DOI:
integer,parameter   ::   XDM_RC_FROM_CX  = 1

! DOI:
integer,parameter   ::   XDM_RC_FROM_POL = 2

integer :: xdm_rc_source    = XDM_RC_FROM_CX

! ------------------------------------------------------------------------------

! https://en.wikipedia.org/wiki/Ionization_energies_of_the_elements_(data_page)
! CRC data in eV

real(DEVDP)     :: xdm_ip(10) = (/  &
13.59844, 24.58741, &
5.39172, 9.3227, 8.29803, 11.26030, 14.53414, 13.61806, 17.42282, 21.5646 &
/)

! ------------------------------------------------------------------------------

end module ffdev_xdm_dat
