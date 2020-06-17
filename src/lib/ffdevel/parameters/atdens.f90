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

module ffdev_atdens

use ffdev_atdens_dat
use ffdev_constants
use ffdev_variables

contains

! ==============================================================================
! function ffdev_atdens_b
! ==============================================================================

real(DEVDP) function ffdev_atdens_b(z)

    use ffdev_utils

    implicit none
    integer     :: z
    ! --------------------------------------------------------------------------

    if( (z .le. 0) .and. (z .gt. ATDENS_MAX_Z) ) then
        call ffdev_utils_exit(DEV_ERR,1,'Z is out-of-range in ffdev_atdens_b')
    end if

    select case(atdens_source)
        case(ATDENS_HF_UGBS)
            if( .not. atdens_HF_UGBS_avail(z) ) then
                call ffdev_utils_exit(DEV_ERR,1,'No data for given Z')
            end if
            ffdev_atdens_b = atdens_HF_UGBS_b(z)

        case(ATDENS_HF_DKH_ANORCC)
            if( .not. atdens_HF_DKH_ANORCC_avail(z) ) then
                call ffdev_utils_exit(DEV_ERR,1,'No data for given Z')
            end if
            ffdev_atdens_b = atdens_HF_DKH_ANORCC_b(z)

        case default
            call ffdev_utils_exit(DEV_ERR,1,'Not implemented in ffdev_atdens_b')
    end select

end function

! ==============================================================================
! function ffdev_atdens_rc
! ==============================================================================

real(DEVDP) function ffdev_atdens_rc(z,dens)

    use ffdev_utils

    implicit none
    integer     :: z
    real(DEVDP) :: dens
    ! --------------------------------------------------------------------------

    if( (z .le. 0) .and. (z .gt. ATDENS_MAX_Z) ) then
        call ffdev_utils_exit(DEV_ERR,1,'Z is out-of-range in ffdev_atdens_rc')
    end if

    select case(atdens_source)
        case(ATDENS_HF_UGBS)
            if( .not. atdens_HF_UGBS_avail(z) ) then
                call ffdev_utils_exit(DEV_ERR,1,'No data for given Z')
            end if
            ffdev_atdens_rc = (atdens_HF_UGBS_a(z) - dens)/atdens_HF_UGBS_b(z)

        case(ATDENS_HF_DKH_ANORCC)
            if( .not. atdens_HF_DKH_ANORCC_avail(z) ) then
                call ffdev_utils_exit(DEV_ERR,1,'No data for given Z')
            end if
            ffdev_atdens_rc = (atdens_HF_DKH_ANORCC_a(z) - dens)/atdens_HF_DKH_ANORCC_b(z)

        case default
            call ffdev_utils_exit(DEV_ERR,1,'Not implemented in ffdev_atdens_rc')
    end select

end function ffdev_atdens_rc

! ==============================================================================
! subroutine ffdev_atdens_source_from_string
! ==============================================================================

integer function ffdev_atdens_source_from_string(string)

    use ffdev_utils

    implicit none
    character(*)   :: string
    ! --------------------------------------------------------------------------

    select case(trim(string))
        case('HF/UGBS')
            ffdev_atdens_source_from_string = ATDENS_HF_UGBS
        case('CC/UGBS')
            ffdev_atdens_source_from_string = ATDENS_CC_UGBS
        case('HF-DKH/ANORCC')
            ffdev_atdens_source_from_string = ATDENS_HF_DKH_ANORCC
        case('CC-DKH/ANORCC')
            ffdev_atdens_source_from_string = ATDENS_CC_DKH_ANORCC
        case default
            call ffdev_utils_exit(DEV_ERR,1,'Not implemented "' // trim(string) //'" in ffdev_atdens_source_from_string!')
    end select

end function ffdev_atdens_source_from_string

! ==============================================================================
! subroutine ffdev_atdens_source_to_string
! ==============================================================================

character(80) function ffdev_atdens_source_to_string(mode)

    use ffdev_utils

    implicit none
    integer  :: mode
    ! --------------------------------------------------------------------------

    select case(mode)
        case(ATDENS_HF_UGBS)
            ffdev_atdens_source_to_string = 'HF/UGBS'
        case(ATDENS_CC_UGBS)
            ffdev_atdens_source_to_string = 'CC/UGBS'
        case(ATDENS_HF_DKH_ANORCC)
            ffdev_atdens_source_to_string = 'HF-DKH/ANORCC'
        case(ATDENS_CC_DKH_ANORCC)
            ffdev_atdens_source_to_string = 'CC-DKH/ANORCC'
        case default
            call ffdev_utils_exit(DEV_ERR,1,'Not implemented in ffdev_atdens_source_to_string!')
    end select

end function ffdev_atdens_source_to_string

! ------------------------------------------------------------------------------

end module ffdev_atdens
