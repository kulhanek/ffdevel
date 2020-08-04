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
        case(DAMP_BJ_RDO)
            ffdev_topology_BJ_dampbj_mode_to_string = 'RDO - Derived from electron density overlaps'
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
        case('RDO')
            ffdev_topology_BJ_dampbj_mode_from_string = DAMP_BJ_RDO
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
    integer         :: i,ti,agti,tj,agtj
    ! --------------------------------------------------------------------------

    ! first update RC from external sources
    select case(dampbj_mode)
        case(DAMP_BJ_CONST)
            do i=1,top%nnb_types
                top%nb_types(i)%tb = damp_fa
            end do
    !---------------
        case(DAMP_BJ_FREEOPT)
            ! nothing
    !---------------
        case(DAMP_BJ_DRC)
            do i=1,top%nnb_types
                ti   = top%nb_types(i)%ti
                tj   = top%nb_types(i)%tj
                agti = top%atom_types(ti)%glbtypeid
                agtj = top%atom_types(tj)%glbtypeid
                top%nb_types(i)%rc = damp_fa * disp_pairs(agti,agtj)%rc + damp_fb
            end do
    !---------------
        case(DAMP_BJ_RDO)
            do i=1,top%nnb_types
                if( top%nb_types(i)%ti .ne. top%nb_types(i)%tj ) cycle
                ti   = top%nb_types(i)%ti
                agti = top%atom_types(ti)%glbtypeid
                top%nb_types(i)%rc = ffdev_atomicdata_rcii(agti,damp_fa)
            end do
    !---------------
        case default
            call ffdev_utils_exit(DEV_ERR,1,'RC mode not implemented in ffdev_topology_BJ_update_nb_params I!')
    end select

    ! apply combining rules if necessary
    select case(dampbj_mode)
        case(DAMP_BJ_CONST,DAMP_BJ_DRC)
            ! nothing
    !---------------
        case(DAMP_BJ_FREEOPT)
            if( ApplyCombiningRules ) then
                call ffdev_topology_BJ_update_nb_params_RC(top)
            end if
    !---------------
        case(DAMP_BJ_RDO)
            if( .not. ApplyCombiningRules ) then
                ! we need to apply combining rules for unlike atoms
                call ffdev_utils_exit(DEV_ERR,1,'DAMP_BJ_DO requires ApplyCombiningRules!')
            end if
            call ffdev_topology_BJ_update_nb_params_RC(top)
    !---------------
        case default
            call ffdev_utils_exit(DEV_ERR,1,'RC mode not implemented in ffdev_topology_BJ_update_nb_params II!')
    end select

end subroutine ffdev_topology_BJ_update_nb_params

!===============================================================================
! subroutine ffdev_topology_BJ_update_nb_params_RC
!===============================================================================

subroutine ffdev_topology_BJ_update_nb_params_RC(top)

    use ffdev_utils
    use ffdev_topology_exp

    implicit none
    type(TOPOLOGY)  :: top
    ! --------------------------------------------
    integer         :: i,nbii,nbjj
    real(DEVDP)     :: rcii,rcij,rcjj
    ! --------------------------------------------------------------------------
    ! apply combining rules - Rc average, only FREEOPT
    do i=1,top%nnb_types

        ! discard like atoms
        if( top%nb_types(i)%ti .eq. top%nb_types(i)%tj ) cycle

        ! get type parameters
        nbii = top%nb_types(i)%nbii
        nbjj = top%nb_types(i)%nbjj

        rcii = top%nb_types(nbii)%rc
        rcjj = top%nb_types(nbjj)%rc

        rcij = 0.5d0 * (rcii+rcjj)

        top%nb_types(i)%rc = rcij
    end do

end subroutine ffdev_topology_BJ_update_nb_params_RC

! ------------------------------------------------------------------------------

end module ffdev_topology_bj

