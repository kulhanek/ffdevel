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


end module ffdev_disp_dat

