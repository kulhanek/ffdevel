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
! bii source
integer,parameter   :: AD_BII_IP                = 1
integer,parameter   :: AD_BII_PBE0_def2QZVPP    = 2

! bii modificators
integer,parameter   :: AD_BII_RAW               = 1     ! raw data
integer,parameter   :: AD_BII_MOD_BY_XDM        = 2     ! modified by XDM volumes
integer,parameter   :: AD_BII_MOD_BY_CHRG       = 3     ! modified by charge

! setup bii
integer             :: bii_source   = AD_BII_IP
integer             :: bii_mods     = AD_BII_RAW

! ------------------------------------------------------------------------------
! ------------------------------------------------------------------------------

! rcii source
integer,parameter   :: AD_RCII_VDW              = 1
integer,parameter   :: AD_RCII_XDM              = 2
integer,parameter   :: AD_RCII_PBE0_def2QZVPP   = 3

! setup bii
integer             :: rcii_source   = AD_RCII_XDM

! ------------------------------------------------------------------------------

end module ffdev_atomicdata_dat
