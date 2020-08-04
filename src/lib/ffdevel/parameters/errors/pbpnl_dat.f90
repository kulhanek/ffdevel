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

module ffdev_err_pbpnl_dat

use ffdev_constants
use ffdev_variables

! ------------------------------------------------------------------------------

! mode
integer,parameter   :: PBPNL_MODE_ALL      = 1
integer,parameter   :: PBPNL_MODE_BURIED   = 2
integer,parameter   :: PBPNL_MODE_NOH      = 3

! PB source
integer,parameter   :: PBPNL_SOURCE_DO     = 1     ! Density overlap
integer,parameter   :: PBPNL_SOURCE_WO     = 2     ! Wavefunction overlap
integer,parameter   :: PBPNL_SOURCE_IP     = 3     ! Ionization potential


! restraint form
integer,parameter   :: PBPNL_QUADRATIC     = 1

! ------------------------------------------------------------------------------

! initialization in ffdev_err_pbpnl_init
logical     :: EnablePBPnlError
logical     :: PrintPBPnlErrorSummary
integer     :: PBPNLMode
integer     :: PBPNLSource
logical     :: PBPNLIncludeProbes

real(DEVDP) :: PBPnlErrorWeight     ! global
real(DEVDP) :: PBPnlErrorWeight1    ! high
real(DEVDP) :: PBPnlErrorWeight2    ! low

integer     :: PBPnlFce

! ------------------------------------------------------------------------------

end module ffdev_err_pbpnl_dat
