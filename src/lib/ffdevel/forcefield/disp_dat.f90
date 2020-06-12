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

module ffdev_disp_dat

use ffdev_sizes
use ffdev_constants
use ffdev_variables

! ==============================================================================
! ==== dispersion
! ==============================================================================

integer,parameter   :: NB_CX_NONE       = 0
integer,parameter   :: NB_CX_XDM        = 1
integer,parameter   :: NB_CX_MMD3       = 2

integer,parameter   :: NB_RC_NONE       = 0
integer,parameter   :: NB_RC_XDM        = 1
integer,parameter   :: NB_RC_XDM_POL    = 2
integer,parameter   :: NB_RC_XDM_VOL    = 3
integer,parameter   :: NB_RC_MMD3       = 4

! dispersion parameters pair
type DISP_PAIR_TYPE
    real(DEVDP)         :: c6
    real(DEVDP)         :: c8
    real(DEVDP)         :: c10
    real(DEVDP)         :: Rc
end type DISP_PAIR_TYPE

type(DISP_PAIR_TYPE),allocatable    :: disp_pairs(:,:)      ! ntypes x ntypes - global types
logical                             :: disp_data_loaded     = .false.
integer                             :: cx_source            = NB_CX_NONE
integer                             :: rc_source            = NB_RC_NONE

! ------------------------------------------------------------------------------

! https://en.wikipedia.org/wiki/Ionization_energies_of_the_elements_(data_page)
! CRC data in eV

real(DEVDP)     :: atom_ip(10) = (/  &
13.59844, 24.58741, &
5.39172, 9.3227, 8.29803, 11.26030, 14.53414, 13.61806, 17.42282, 21.5646 &
/)


! ------------------------------------------------------------------------------

! CCSD/UKBS - bfac

real(DEVDP)     :: atom_bfac(18) = (/  &
3.7829, 5.3139, &
1.9827, 2.9793, 3.6458, 3.4357, 4.2317, 4.7252, 4.6843, 5.2192, &
1.8440, 2.5588, 2.3121, 2.7842, 3.2578, 3.2282, 4.2113, 4.3007 &
/)


end module ffdev_disp_dat

