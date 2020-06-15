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

module ffdev_disp

use ffdev_sizes
use ffdev_constants
use ffdev_variables

contains

! ==============================================================================
! subroutine ffdev_disp_cxsource_from_string
! ==============================================================================

integer function ffdev_disp_cxsource_from_string(string)

    use ffdev_utils
    use ffdev_disp_dat

    implicit none
    character(*)   :: string
    ! --------------------------------------------------------------------------

    select case(trim(string))
        case('NONE')
            ffdev_disp_cxsource_from_string = NB_CX_NONE
        case('XDM')
            ffdev_disp_cxsource_from_string = NB_CX_XDM
        case('MMD3')
            ffdev_disp_cxsource_from_string = NB_CX_MMD3
        case default
            call ffdev_utils_exit(DEV_ERR,1,'Not implemented "' // trim(string) //'" in ffdev_disp_cxsource_from_string!')
    end select

end function ffdev_disp_cxsource_from_string

! ==============================================================================
! subroutine ffdev_disp_cxsource_to_string
! ==============================================================================

character(80) function ffdev_disp_cxsource_to_string(nb_mode)

    use ffdev_utils
    use ffdev_disp_dat

    implicit none
    integer  :: nb_mode
    ! --------------------------------------------------------------------------

    select case(nb_mode)
        case(NB_CX_NONE)
            ffdev_disp_cxsource_to_string = 'NONE'
        case(NB_CX_XDM)
            ffdev_disp_cxsource_to_string = 'XDM'
        case(NB_CX_MMD3)
            ffdev_disp_cxsource_to_string = 'MMD3'
        case default
            call ffdev_utils_exit(DEV_ERR,1,'Not implemented in ffdev_disp_cxsource_to_string!')
    end select

end function ffdev_disp_cxsource_to_string

! ==============================================================================
! subroutine ffdev_disp_rcsource_from_string
! ==============================================================================

integer function ffdev_disp_rcsource_from_string(string)

    use ffdev_utils
    use ffdev_disp_dat

    implicit none
    character(*)   :: string
    ! --------------------------------------------------------------------------

    select case(trim(string))
        case('NONE')
            ffdev_disp_rcsource_from_string = NB_RC_NONE
        case('XDM')
            ffdev_disp_rcsource_from_string = NB_RC_XDM
        case('XDM-POL')
            ffdev_disp_rcsource_from_string = NB_RC_XDM_POL
        case('XDM-VOL')
            ffdev_disp_rcsource_from_string = NB_RC_XDM_VOL
        case('MMD3')
            ffdev_disp_rcsource_from_string = NB_RC_MMD3
        case default
            call ffdev_utils_exit(DEV_ERR,1,'Not implemented "' // trim(string) //'" in ffdev_disp_rcsource_from_string!')
    end select

end function ffdev_disp_rcsource_from_string

! ==============================================================================
! subroutine ffdev_disp_rcsource_to_string
! ==============================================================================

character(80) function ffdev_disp_rcsource_to_string(nb_mode)

    use ffdev_utils
    use ffdev_disp_dat

    implicit none
    integer  :: nb_mode
    ! --------------------------------------------------------------------------

    select case(nb_mode)
        case(NB_RC_NONE)
            ffdev_disp_rcsource_to_string = 'NONE'
        case(NB_RC_XDM)
            ffdev_disp_rcsource_to_string = 'XDM'
        case(NB_RC_XDM_POL)
            ffdev_disp_rcsource_to_string = 'XDM-POL'
        case(NB_RC_XDM_VOL)
            ffdev_disp_rcsource_to_string = 'XDM-VOL'
        case(NB_RC_MMD3)
            ffdev_disp_rcsource_to_string = 'MMD3'
        case default
            call ffdev_utils_exit(DEV_ERR,1,'Not implemented in ffdev_disp_rcsource_to_string!')
    end select

end function ffdev_disp_rcsource_to_string

! ==============================================================================

end module ffdev_disp

