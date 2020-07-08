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

module ffdev_errors_utils

use ffdev_constants
use ffdev_variables

contains

! ==============================================================================
! subroutine ffdev_errors_utils_scale_to_string
! ==============================================================================

character(MAX_PATH) function ffdev_errors_utils_scale_to_string(iscale)

    use ffdev_utils
    use ffdev_errors_dat

    implicit none
    integer  :: iscale
    ! --------------------------------------------------------------------------

    select case(iscale)
        case(EE_ABS)
            ffdev_errors_utils_scale_to_string = 'abs (absolute)'
        case(EE_REL)
            ffdev_errors_utils_scale_to_string = 'rel (relative)'
        case(EE_LOG)
            ffdev_errors_utils_scale_to_string = 'log (logarithimic)'
        case default
            call ffdev_utils_exit(DEV_ERR,1,'Not implemented in ffdev_errors_utils_scale_to_string!')
    end select

end function ffdev_errors_utils_scale_to_string

! ==============================================================================
! function ffdev_errors_utils_scale_from_string
! ==============================================================================

integer function ffdev_errors_utils_scale_from_string(string)

    use ffdev_utils
    use ffdev_errors_dat

    implicit none
    character(*)   :: string
    ! --------------------------------------------------------------------------

    select case(trim(string))
        case('abs')
            ffdev_errors_utils_scale_from_string = EE_ABS
        case('rel')
            ffdev_errors_utils_scale_from_string = EE_REL
        case('log')
            ffdev_errors_utils_scale_from_string = EE_LOG
        case default
            call ffdev_utils_exit(DEV_ERR,1,'Not implemented "' // trim(string) // &
                                            '" in ffdev_errors_utils_scale_from_string!')
    end select

end function ffdev_errors_utils_scale_from_string

! ------------------------------------------------------------------------------

end module ffdev_errors_utils
