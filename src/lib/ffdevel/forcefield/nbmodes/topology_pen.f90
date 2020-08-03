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

module ffdev_topology_pen

use ffdev_topology_dat
use ffdev_constants
use ffdev_variables

contains

! ==============================================================================
! subroutine ffdev_topology_PEN_pen_mode_to_string
! ==============================================================================

character(80) function ffdev_topology_PEN_pen_mode_to_string(lexp_mode)

    use ffdev_utils

    implicit none
    integer  :: lexp_mode
    ! --------------------------------------------------------------------------

    select case(lexp_mode)
        case(PEN_MODE1)
            ffdev_topology_PEN_pen_mode_to_string = 'MODE1 - test1'
        case default
            call ffdev_utils_exit(DEV_ERR,1,'Not implemented in ffdev_topology_PEN_pen_mode_to_string!')
    end select

end function ffdev_topology_PEN_pen_mode_to_string

! ==============================================================================
! subroutine ffdev_topology_PEN_pen_mode_from_string
! ==============================================================================

integer function ffdev_topology_PEN_pen_mode_from_string(string)

    use ffdev_utils

    implicit none
    character(*)   :: string
    ! --------------------------------------------------------------------------

    select case(trim(string))
        case('MODE1')
            ffdev_topology_PEN_pen_mode_from_string = PEN_MODE1
        case default
            call ffdev_utils_exit(DEV_ERR,1,'Not implemented "' // trim(string) //'" in ffdev_topology_PEN_pen_mode_from_string!')
    end select

end function ffdev_topology_PEN_pen_mode_from_string

! ------------------------------------------------------------------------------

end module ffdev_topology_pen

