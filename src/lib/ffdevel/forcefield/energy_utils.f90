! ==============================================================================
! This file is part of FFDevel.
!    Copyright (C) 2018 Petr Kulhanek, kulhanek@chemi.muni.cz
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

module ffdev_energy_utils

use ffdev_geometry_dat
use ffdev_constants

contains

! ==============================================================================
! subroutine ffdev_energy_utils_allocate_eneprmgrad
! ==============================================================================

subroutine ffdev_energy_utils_allocate_eneprmgrad(geo)

    use ffdev_geometry
    use ffdev_utils
    use ffdev_parameters_dat

    implicit none
    type(GEOMETRY)  :: geo
    ! --------------------------------------------
    integer         :: alloc_stat
    ! --------------------------------------------------------------------------

    if( nparams .le. 0 ) then
        call ffdev_utils_exit(DEV_OUT,1,'Number of parameters has to be greater than zero for eneprmgrd!')
    end if

    allocate(geo%eneprmgrd(nparams), stat = alloc_stat)
    if( alloc_stat .ne. 0 ) then
        call ffdev_utils_exit(DEV_OUT,1,'Unable to allocate temporary array for eneprmgrd!')
    end if

    geo%eneprmgrd(:) = 0.0d0

end subroutine ffdev_energy_utils_allocate_eneprmgrad

! ------------------------------------------------------------------------------

end module ffdev_energy_utils
