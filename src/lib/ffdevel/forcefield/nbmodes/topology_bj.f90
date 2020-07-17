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
        case(DAMP_BJ_DRC)
            ffdev_topology_BJ_dampbj_mode_to_string = 'DRC - Use Rc from dispersion data'
        case(DAMP_BJ_DO)
            ffdev_topology_BJ_dampbj_mode_to_string = 'DO - Derived from electron density overlaps'
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
        case('DRC')
            ffdev_topology_BJ_dampbj_mode_from_string = DAMP_BJ_DRC
        case('DO')
            ffdev_topology_BJ_dampbj_mode_from_string = DAMP_BJ_DO
        case default
            call ffdev_utils_exit(DEV_ERR,1,'Not implemented "' // trim(string) //'" in ffdev_topology_BJ_dampbj_mode_from_string!')
    end select

end function ffdev_topology_BJ_dampbj_mode_from_string

! ==============================================================================
! subroutine ffdev_topology_BJ_apply_NB_comb_rules
! ==============================================================================

subroutine ffdev_topology_BJ_apply_NB_comb_rules(top)

    use ffdev_utils
    use ffdev_topology_utils

    implicit none
    type(TOPOLOGY)  :: top
    ! --------------------------------------------
    integer         :: i,nbii,nbjj
    real(DEVDP)     :: rcii,rcij,rcjj
    ! --------------------------------------------------------------------------

    if( dampbj_mode .ne. DAMP_BJ_FREEOPT ) return

    ! apply combining rules - Rc average, only FREEOPT
    do i=1,top%nnb_types

        ! discard like atoms
        if( top%nb_types(i)%ti .eq. top%nb_types(i)%tj ) cycle

        ! get type parameters
        nbii = ffdev_topology_find_nbtype_by_tindex(top,top%nb_types(i)%ti,top%nb_types(i)%ti)
        nbjj = ffdev_topology_find_nbtype_by_tindex(top,top%nb_types(i)%tj,top%nb_types(i)%tj)

        rcii = top%nb_types(nbii)%rc
        rcjj = top%nb_types(nbjj)%rc

        rcij = 0.5d0 * (rcii+rcjj)

        top%nb_types(i)%rc = rcij
    end do

end subroutine ffdev_topology_BJ_apply_NB_comb_rules

! ------------------------------------------------------------------------------

end module ffdev_topology_bj

