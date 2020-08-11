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
integer,parameter   :: AD_BII_IPEA              = 1     ! ionization potentials / electron affinities
integer,parameter   :: AD_BII_ATDENS            = 2     ! spherically symmetrized atom densities

! bii modificators
integer,parameter   :: AD_BII_RAW               = 1     ! raw data
integer,parameter   :: AD_BII_MOD_BY_XDM        = 2     ! modified by XDM volumes
integer,parameter   :: AD_BII_MOD_BY_CHRG       = 3     ! modified by charge

! atdens source
integer,parameter   :: AD_ATDENS_PBE0_def2QZVPP = 1     ! PBE0/def2-QZVPP

! setup bii
integer             :: bii_source       = AD_BII_IPEA
integer             :: bii_mods         = AD_BII_RAW
integer             :: atdens_source    = AD_ATDENS_PBE0_def2QZVPP

! ------------------------------------------------------------------------------
! ------------------------------------------------------------------------------

! nuclear core

integer,parameter   :: AD_ZEFF_MAX              = 1     ! all electrons
integer,parameter   :: AD_ZEFF_VALENCE          = 2     ! valence electrons
integer,parameter   :: AD_ZEFF_OPT              = 3     ! optimized

integer             :: Zeff_mode                = AD_ZEFF_VALENCE

! ------------------------------------------------------------------------------
! ------------------------------------------------------------------------------

! rcii source
integer,parameter   :: AD_RCII_VDW              = 1
integer,parameter   :: AD_RCII_DISP             = 2
integer,parameter   :: AD_RCII_ATDENS           = 3

! setup bii
integer             :: rcii_source   = AD_RCII_DISP

! ------------------------------------------------------------------------------

end module ffdev_atomicdata_dat
