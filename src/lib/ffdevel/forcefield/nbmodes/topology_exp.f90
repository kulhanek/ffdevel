! ==============================================================================
! This file is part of FFDevel.
!    Copyright (C) 2013 Petr Kulhanek, kulhanek@chemi.muni.cz
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

module ffdev_topology_exp

use ffdev_topology_dat
use ffdev_constants
use ffdev_variables

contains

! ==============================================================================
! subroutine ffdev_topology_EXP_exp_mode_to_string
! ==============================================================================

character(80) function ffdev_topology_EXP_exp_mode_to_string(lexp_mode)

    use ffdev_utils

    implicit none
    integer  :: lexp_mode
    ! --------------------------------------------------------------------------

    select case(lexp_mode)
        case(EXP_MODE_DO)
            ffdev_topology_EXP_exp_mode_to_string = 'EXP-DO - Density overlap'
        case(EXP_MODE_WO)
            ffdev_topology_EXP_exp_mode_to_string = 'EXP-WO - Wavefunction overlap'
        case default
            call ffdev_utils_exit(DEV_ERR,1,'Not implemented in ffdev_topology_EXP_exp_mode_to_string!')
    end select

end function ffdev_topology_EXP_exp_mode_to_string

! ==============================================================================
! subroutine ffdev_topology_EXP_exp_mode_from_string
! ==============================================================================

integer function ffdev_topology_EXP_exp_mode_from_string(string)

    use ffdev_utils

    implicit none
    character(*)   :: string
    ! --------------------------------------------------------------------------

    select case(trim(string))
        case('EXP-DO')
            ffdev_topology_EXP_exp_mode_from_string = EXP_MODE_DO
        case('EXP-WO')
            ffdev_topology_EXP_exp_mode_from_string = EXP_MODE_WO
        case default
            call ffdev_utils_exit(DEV_ERR,1,'Not implemented "' // trim(string) //'" in ffdev_topology_EXP_exp_mode_from_string!')
    end select

end function ffdev_topology_EXP_exp_mode_from_string

! ==============================================================================
! subroutine ffdev_topology_EXP_pa_mode_to_string
! ==============================================================================

character(80) function ffdev_topology_EXP_pa_mode_to_string(lpb_mode)

    use ffdev_utils

    implicit none
    integer  :: lpb_mode
    ! --------------------------------------------------------------------------

    select case(lpb_mode)
        case(EXP_PA_FREEOPT)
            ffdev_topology_EXP_pa_mode_to_string = 'FREEOPT - Use optimized PA per type'
        case(EXP_PA_CHARGES)
            ffdev_topology_EXP_pa_mode_to_string = 'CHARGES - Use effective charges and k_exc'
        case default
            call ffdev_utils_exit(DEV_ERR,1,'Not implemented in ffdev_topology_EXP_pa_mode_to_string!')
    end select

end function ffdev_topology_EXP_pa_mode_to_string

! ==============================================================================
! subroutine ffdev_topology_EXP_pa_mode_from_string
! ==============================================================================

integer function ffdev_topology_EXP_pa_mode_from_string(string)

    use ffdev_utils

    implicit none
    character(*)   :: string
    ! --------------------------------------------------------------------------

    select case(trim(string))
        case('FREEOPT')
            ffdev_topology_EXP_pa_mode_from_string = EXP_PA_FREEOPT
        case('CHARGES')
            ffdev_topology_EXP_pa_mode_from_string = EXP_PA_CHARGES
        case default
            call ffdev_utils_exit(DEV_ERR,1,'Not implemented "' // trim(string) //'" in ffdev_topology_EXP_pa_mode_from_string!')
    end select

end function ffdev_topology_EXP_pa_mode_from_string

! ==============================================================================
! subroutine ffdev_topology_EXP_pb_mode_to_string
! ==============================================================================

character(80) function ffdev_topology_EXP_pb_mode_to_string(lpb_mode)

    use ffdev_utils

    implicit none
    integer  :: lpb_mode
    ! --------------------------------------------------------------------------

    select case(lpb_mode)
        case(EXP_PB_FREEOPT)
            ffdev_topology_EXP_pb_mode_to_string = 'FREEOPT - Use optimized PB per type'
        case(EXP_PB_ADBII)
            ffdev_topology_EXP_pb_mode_to_string = 'ADBII - Atomic database Bii'
        case default
            call ffdev_utils_exit(DEV_ERR,1,'Not implemented in ffdev_topology_EXP_pb_mode_to_string!')
    end select

end function ffdev_topology_EXP_pb_mode_to_string

! ==============================================================================
! subroutine ffdev_topology_EXP_pb_mode_from_string
! ==============================================================================

integer function ffdev_topology_EXP_pb_mode_from_string(string)

    use ffdev_utils

    implicit none
    character(*)   :: string
    ! --------------------------------------------------------------------------

    select case(trim(string))
        case('FREEOPT')
            ffdev_topology_EXP_pb_mode_from_string = EXP_PB_FREEOPT
        case('ADBII')
            ffdev_topology_EXP_pb_mode_from_string = EXP_PB_ADBII
        case default
            call ffdev_utils_exit(DEV_ERR,1,'Not implemented "' // trim(string) //'" in ffdev_topology_EXP_pb_mode_from_string!')
    end select

end function ffdev_topology_EXP_pb_mode_from_string

! ==============================================================================
! subroutine ffdev_topology_EXP_update_nb_params
! ==============================================================================

subroutine ffdev_topology_EXP_update_nb_params(top)

    use ffdev_utils
    use ffdev_atomicdata

    implicit none
    type(TOPOLOGY)  :: top
    ! --------------------------------------------
    integer         :: i,agti
    ! --------------------------------------------------------------------------

    ! first update PB from external sources
    select case(exp_pb_mode)
        case(EXP_PB_FREEOPT)
            ! nothing to do
    !---------------
        case(EXP_PB_ADBII)
            do i=1,top%natom_types
                agti = top%atom_types(i)%glbtypeid
                top%atom_types(i)%pb = ffdev_atomicdata_bii(agti)
            end do
    !---------------
        case default
            call ffdev_utils_exit(DEV_ERR,1,'exp_pb_mode not implemented in ffdev_topology_EXP_update_nb_params I!')
    end select

end subroutine ffdev_topology_EXP_update_nb_params

! ------------------------------------------------------------------------------

end module ffdev_topology_exp

