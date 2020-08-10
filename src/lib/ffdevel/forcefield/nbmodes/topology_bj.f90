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

module ffdev_topology_bj

use ffdev_topology_dat
use ffdev_constants
use ffdev_variables

contains

! ==============================================================================
! subroutine ffdev_topology_dampbj_mode_to_string
! ==============================================================================

character(80) function ffdev_topology_BJ_dampbj_mode_to_string(nb_mode)

    use ffdev_utils

    implicit none
    integer  :: nb_mode
    ! --------------------------------------------------------------------------

    select case(nb_mode)
        case(DAMP_BJ_CONST)
            ffdev_topology_BJ_dampbj_mode_to_string = 'CONST - Constant Rc'
        case(DAMP_BJ_FREEOPT)
            ffdev_topology_BJ_dampbj_mode_to_string = 'FREEOPT - Use optimized Rc per type'
        case(DAMP_BJ_ADRCII)
            ffdev_topology_BJ_dampbj_mode_to_string = 'ADRCII - Atomic database Rcii'
        case default
            call ffdev_utils_exit(DEV_ERR,1,'Not implemented in ffdev_topology_BJ_dampbj_mode_to_string!')
    end select

end function ffdev_topology_BJ_dampbj_mode_to_string

! ==============================================================================
! subroutine ffdev_topology_BJ_dampbj_mode_from_string
! ==============================================================================

integer function ffdev_topology_BJ_dampbj_mode_from_string(string)

    use ffdev_utils

    implicit none
    character(*)   :: string
    ! --------------------------------------------------------------------------

    select case(trim(string))
        case('CONST')
            ffdev_topology_BJ_dampbj_mode_from_string = DAMP_BJ_CONST
        case('FREEOPT')
            ffdev_topology_BJ_dampbj_mode_from_string = DAMP_BJ_FREEOPT
        case('ADRCII')
            ffdev_topology_BJ_dampbj_mode_from_string = DAMP_BJ_ADRCII
        case default
            call ffdev_utils_exit(DEV_ERR,1,'Not implemented "' // trim(string) //'" in ffdev_topology_BJ_dampbj_mode_from_string!')
    end select

end function ffdev_topology_BJ_dampbj_mode_from_string

! ==============================================================================
! subroutine ffdev_topology_BJ_update_nb_params
! ==============================================================================

subroutine ffdev_topology_BJ_update_nb_params(top)

    use ffdev_utils
    use ffdev_disp_dat
    use ffdev_atomicdata

    implicit none
    type(TOPOLOGY)  :: top
    ! --------------------------------------------
    integer         :: i,agti
    ! --------------------------------------------------------------------------

    ! first update RC from external sources
    select case(dampbj_mode)
        case(DAMP_BJ_CONST)
            do i=1,top%natom_types
                agti = top%atom_types(i)%glbtypeid
                top%atom_types(i)%rc = damp_fa
            end do
    !---------------
        case(DAMP_BJ_FREEOPT)
            ! nothing
    !---------------
        case(DAMP_BJ_ADRCII)
            do i=1,top%natom_types
                agti = top%atom_types(i)%glbtypeid
                top%atom_types(i)%rc = ffdev_atomicdata_rcii(agti,damp_fa,damp_fb)
            end do
    !---------------
        case default
            call ffdev_utils_exit(DEV_ERR,1,'RC mode not implemented in ffdev_topology_BJ_update_nb_params I!')
    end select

end subroutine ffdev_topology_BJ_update_nb_params

! ------------------------------------------------------------------------------

end module ffdev_topology_bj

