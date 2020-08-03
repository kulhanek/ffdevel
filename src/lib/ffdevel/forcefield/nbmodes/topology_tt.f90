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
            ffdev_topology_TT_damptt_mode_to_string = 'CONST - a single constant TB'
        case(DAMP_TT_FREEOPT)
            ffdev_topology_TT_damptt_mode_to_string = 'FREEOPT - Optimized TB per type'
        case(DAMP_TT_COUPLED)
            ffdev_topology_TT_damptt_mode_to_string = 'COUPLED - TB coupled by damp_tb to PB'
        case(DAMP_TT_DO)
            ffdev_topology_TT_damptt_mode_to_string = 'DO - TB derived from density overlap'
        case(DAMP_TT_IP)
            ffdev_topology_TT_damptt_mode_to_string = 'IP - TB derived from ionization potentials'
        case(DAMP_TT_IP_XDM)
            ffdev_topology_TT_damptt_mode_to_string = 'IP-XDM - TB derived from ionization potentials + XDM mods'
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
        case('IP')
            ffdev_topology_TT_damptt_mode_from_string = DAMP_TT_IP
        case('IP-XDM')
            ffdev_topology_TT_damptt_mode_from_string = DAMP_TT_IP_XDM
        case default
            call ffdev_utils_exit(DEV_ERR,1,'Not implemented "' // trim(string) //'" in ffdev_topology_TT_damptt_mode_from_string!')
    end select

end function ffdev_topology_TT_damptt_mode_from_string

! ==============================================================================
! subroutine ffdev_topology_TT_update_nb_params
! ==============================================================================

subroutine ffdev_topology_TT_update_nb_params(top)

    use ffdev_utils
    use ffdev_atomoverlap
    use ffdev_ip_db

    implicit none
    type(TOPOLOGY)  :: top
    ! --------------------------------------------
    integer         :: i,ti,agti,tj,agtj
    ! --------------------------------------------------------------------------

    ! first update TB from external sources
    select case(damptt_mode)
        case(DAMP_TT_CONST)
            do i=1,top%nnb_types
                top%nb_types(i)%tb = 1.0d0
            end do
    !---------------
        case(DAMP_TT_FREEOPT)
            ! nothing
    !---------------
        case(DAMP_TT_COUPLED)
            do i=1,top%nnb_types
                top%nb_types(i)%tb = top%nb_types(i)%pb
            end do
    !---------------
        case(DAMP_TT_DO)
            do i=1,top%nnb_types
                if( top%nb_types(i)%ti .ne. top%nb_types(i)%tj ) cycle
                ti   = top%nb_types(i)%ti
                agti = top%atom_types(ti)%glbtypeid
                top%nb_types(i)%tb = ffdev_atomoverlap_do_bii(agti)
            end do
    !---------------
        case(DAMP_TT_IP)
            do i=1,top%nnb_types
                if( top%nb_types(i)%ti .ne. top%nb_types(i)%tj ) cycle
                ti   = top%nb_types(i)%ti
                agti = top%atom_types(ti)%glbtypeid
                top%nb_types(i)%tb = ffdev_bfac_from_ip(agti)
            end do
    !---------------
        case(DAMP_TT_IP_XDM)
            do i=1,top%nnb_types
                if( top%nb_types(i)%ti .ne. top%nb_types(i)%tj ) cycle
                ti   = top%nb_types(i)%ti
                agti = top%atom_types(ti)%glbtypeid
                top%nb_types(i)%tb = ffdev_bfac_from_ip_xdm(agti)
            end do
    !---------------
        case default
            call ffdev_utils_exit(DEV_ERR,1,'TT damp mode not implemented in ffdev_topology_TT_update_nb_params I!')
    end select

    ! apply combining rules if necessary
    select case(damptt_mode)
        case(DAMP_TT_CONST,DAMP_TT_COUPLED)
            ! nothing
    !---------------
        case(DAMP_TT_FREEOPT)
            if( ApplyCombiningRules ) then
                call ffdev_topology_TT_update_nb_params_TB(top)
            end if
    !---------------
        case(DAMP_TT_DO,DAMP_TT_IP,DAMP_TT_IP_XDM)
            if( .not. ApplyCombiningRules ) then
                ! we need to apply combining rules for unlike atoms
                call ffdev_utils_exit(DEV_ERR,1,'DAMP_TT_DO, DAMP_TT_IP, or DAMP_TT_IP_XDM requires ApplyCombiningRules!')
            end if
            call ffdev_topology_TT_update_nb_params_TB(top)
    !---------------
        case default
            call ffdev_utils_exit(DEV_ERR,1,'TT damp mode not implemented in ffdev_topology_TT_update_nb_params II!')
    end select

end subroutine ffdev_topology_TT_update_nb_params

!===============================================================================
! subroutine ffdev_topology_TT_update_nb_params_TB
!===============================================================================

subroutine ffdev_topology_TT_update_nb_params_TB(top)

    use ffdev_utils
    use ffdev_topology_exp

    implicit none
    type(TOPOLOGY)  :: top
    ! --------------------------------------------
    integer         :: i,nbii,nbjj
    real(DEVDP)     :: tbii,tbij,tbjj
    ! --------------------------------------------------------------------------

    ! apply combining rules, only FREEOPT
    do i=1,top%nnb_types

        ! discard like atoms
        if( top%nb_types(i)%ti .eq. top%nb_types(i)%tj ) cycle

        ! get type parameters
        nbii = top%nb_types(i)%nbii
        nbjj = top%nb_types(i)%nbjj

        tbii = top%nb_types(nbii)%tb
        tbjj = top%nb_types(nbjj)%tb

        call ffdev_topology_EXP_apply_NB_comb_rules_PB(tbii,tbjj,tbij)

        top%nb_types(i)%tb = tbij
    end do

end subroutine ffdev_topology_TT_update_nb_params_TB

! ------------------------------------------------------------------------------

end module ffdev_topology_tt

