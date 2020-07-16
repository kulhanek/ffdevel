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

module ffdev_topology_tt

use ffdev_topology_dat
use ffdev_constants
use ffdev_variables

contains

! ==============================================================================
! subroutine ffdev_topology_TT_damptt_mode_to_string
! ==============================================================================

character(80) function ffdev_topology_TT_damptt_mode_to_string(nb_mode)

    use ffdev_utils

    implicit none
    integer  :: nb_mode
    ! --------------------------------------------------------------------------

    select case(nb_mode)
        case(DAMP_TT_CONST)
            ffdev_topology_TT_damptt_mode_to_string = 'CONST - constant'
        case(DAMP_TT_FREEOPT)
            ffdev_topology_TT_damptt_mode_to_string = 'FREEOPT - Use optimized TB per type'
        case(DAMP_TT_COUPLED)
            ffdev_topology_TT_damptt_mode_to_string = 'COUPLED - TB coupled by damp_pb to PB'
        case(DAMP_TT_DO)
            ffdev_topology_TT_damptt_mode_to_string = 'DO - b-factors form density overlap'
        case default
            call ffdev_utils_exit(DEV_ERR,1,'Not implemented in ffdev_topology_TT_damptt_mode_to_string!')
    end select

end function ffdev_topology_TT_damptt_mode_to_string

! ==============================================================================
! subroutine ffdev_topology_TT_damptt_mode_from_string
! ==============================================================================

integer function ffdev_topology_TT_damptt_mode_from_string(string)

    use ffdev_utils
    use ffdev_xdm_dat

    implicit none
    character(*)   :: string
    ! --------------------------------------------------------------------------

    select case(trim(string))
        case('CONST')
            ffdev_topology_TT_damptt_mode_from_string = DAMP_TT_CONST
        case('FREEOPT')
            ffdev_topology_TT_damptt_mode_from_string = DAMP_TT_FREEOPT
        case('COUPLED')
            ffdev_topology_TT_damptt_mode_from_string = DAMP_TT_COUPLED
        case('DO')
            ffdev_topology_TT_damptt_mode_from_string = DAMP_TT_DO
        case default
            call ffdev_utils_exit(DEV_ERR,1,'Not implemented "' // trim(string) //'" in ffdev_topology_TT_damptt_mode_from_string!')
    end select

end function ffdev_topology_TT_damptt_mode_from_string

! ==============================================================================
! subroutine ffdev_topology_TT_apply_NB_comb_rules
! ==============================================================================

subroutine ffdev_topology_TT_apply_NB_comb_rules(top)

    use ffdev_utils
    use ffdev_topology_utils
    use ffdev_topology_exp

    implicit none
    type(TOPOLOGY)  :: top
    ! --------------------------------------------
    integer         :: i,nbii,nbjj
    real(DEVDP)     :: tbii,tbij,tbjj
    ! --------------------------------------------------------------------------

    if( damptt_mode .ne. DAMP_TT_FREEOPT ) return

    ! apply combining rules, only FREEOPT
    do i=1,top%nnb_types
        if( top%nb_types(i)%ti .ne. top%nb_types(i)%tj ) then

            ! get type parameters
            nbii = ffdev_topology_find_nbtype_by_tindex(top,top%nb_types(i)%ti,top%nb_types(i)%ti)
            nbjj = ffdev_topology_find_nbtype_by_tindex(top,top%nb_types(i)%tj,top%nb_types(i)%tj)

            tbii = top%nb_types(nbii)%tb
            tbjj = top%nb_types(nbjj)%tb

            call ffdev_topology_EXP_apply_NB_comb_rules_PB(tbii,tbjj,tbij)

            top%nb_types(i)%tb = tbij

        end if
    end do

end subroutine ffdev_topology_TT_apply_NB_comb_rules

! ------------------------------------------------------------------------------

end module ffdev_topology_tt

