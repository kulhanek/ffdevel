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

module ffdev_err_aimr0_dat

use ffdev_constants
use ffdev_variables

! ------------------------------------------------------------------------------

! initialization in ffdev_err_aimr0_init
logical     :: EnableAIMR0Error
logical     :: PrintAIMR0ErrorSummary
real(DEVDP) :: AIMR0ErrorWeight
logical     :: AIMR0BurriedOnly

! ------------------------------------------------------------------------------

end module ffdev_err_aimr0_dat
