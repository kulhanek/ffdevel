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
        case(EXP_BM)
            ffdev_topology_EXP_exp_mode_to_string = 'EXP-BM - Born-Mayer model'
        case(EXP_DO)
            ffdev_topology_EXP_exp_mode_to_string = 'EXP-DO - Density overlap'
        case(EXP_WO)
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
        case('EXP-BM')
            ffdev_topology_EXP_exp_mode_from_string = EXP_BM
        case('EXP-DO')
            ffdev_topology_EXP_exp_mode_from_string = EXP_DO
        case('EXP-WO')
            ffdev_topology_EXP_exp_mode_from_string = EXP_WO
        case default
            call ffdev_utils_exit(DEV_ERR,1,'Not implemented "' // trim(string) //'" in ffdev_topology_EXP_exp_mode_from_string!')
    end select

end function ffdev_topology_EXP_exp_mode_from_string

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
        case(EXP_PB_DO)
            ffdev_topology_EXP_pb_mode_to_string = 'DO - PB from density overlap'
        case(EXP_PB_WO)
            ffdev_topology_EXP_pb_mode_to_string = 'WO - PB from wavefunction overlap'
        case(EXP_PB_IP)
            ffdev_topology_EXP_pb_mode_to_string = 'IP - PB derived from ionization potentials'
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
        case('DO')
            ffdev_topology_EXP_pb_mode_from_string = EXP_PB_DO
        case('WO')
            ffdev_topology_EXP_pb_mode_from_string = EXP_PB_WO
        case('IP')
            ffdev_topology_EXP_pb_mode_from_string = EXP_PB_IP
        case default
            call ffdev_utils_exit(DEV_ERR,1,'Not implemented "' // trim(string) //'" in ffdev_topology_EXP_pb_mode_from_string!')
    end select

end function ffdev_topology_EXP_pb_mode_from_string

! ==============================================================================
! subroutine ffdev_topology_EXP_comb_rules_to_string
! ==============================================================================

character(80) function ffdev_topology_EXP_comb_rules_to_string(comb_rules)

    use ffdev_utils

    implicit none
    integer  :: comb_rules
    ! --------------------------------------------------------------------------

    select case(comb_rules)

        ! EXP potential
        case(EXP_COMB_RULE_AM)
            ffdev_topology_EXP_comb_rules_to_string = 'AM (arithmetic means)'
        case(EXP_COMB_RULE_GS)
            ffdev_topology_EXP_comb_rules_to_string = 'GS (Gilbert-Smith)'
        case(EXP_COMB_RULE_BA)
            ffdev_topology_EXP_comb_rules_to_string = 'BA (Bohm-Ahlrichs)'
        case(EXP_COMB_RULE_VS)
            ffdev_topology_EXP_comb_rules_to_string = 'VS (Vleet-Schmidt)'
        case(EXP_COMB_RULE_D1)
            ffdev_topology_EXP_comb_rules_to_string = 'D1 (PB - Density Overlap, PA=AM)'
        case(EXP_COMB_RULE_D2)
            ffdev_topology_EXP_comb_rules_to_string = 'D2 (PA+PB - Density Overlap)'
        case(EXP_COMB_RULE_W1)
            ffdev_topology_EXP_comb_rules_to_string = 'W1 (PB - Wavefunction Overlap, PA=AM)'
        case(EXP_COMB_RULE_W2)
            ffdev_topology_EXP_comb_rules_to_string = 'W2 (PA+PB - Wavefunction Overlap)'
        case default
            call ffdev_utils_exit(DEV_ERR,1,'Not implemented in ffdev_topology_EXP_comb_rules_to_string!')
    end select

end function ffdev_topology_EXP_comb_rules_to_string

! ==============================================================================
! function ffdev_topology_EXP_comb_rules_from_string
! ==============================================================================

integer function ffdev_topology_EXP_comb_rules_from_string(string)

    use ffdev_utils

    implicit none
    character(*)   :: string
    ! --------------------------------------------------------------------------

    select case(trim(string))
        case('AM')
            ffdev_topology_EXP_comb_rules_from_string = EXP_COMB_RULE_AM
        case('GS')
            ffdev_topology_EXP_comb_rules_from_string = EXP_COMB_RULE_GS
        case('BA')
            ffdev_topology_EXP_comb_rules_from_string = EXP_COMB_RULE_BA
        case('VS')
            ffdev_topology_EXP_comb_rules_from_string = EXP_COMB_RULE_VS
        case('D1')
            ffdev_topology_EXP_comb_rules_from_string = EXP_COMB_RULE_D1
        case('D2')
            ffdev_topology_EXP_comb_rules_from_string = EXP_COMB_RULE_D2
        case('W1')
            ffdev_topology_EXP_comb_rules_from_string = EXP_COMB_RULE_W1
        case('W2')
            ffdev_topology_EXP_comb_rules_from_string = EXP_COMB_RULE_W2
        case default
            call ffdev_utils_exit(DEV_ERR,1,'Not implemented "' // trim(string) //'" in ffdev_topology_EXP_comb_rules_from_string!')
    end select

end function ffdev_topology_EXP_comb_rules_from_string

! ==============================================================================
! subroutine ffdev_topology_EXP_update_nb_params
! ==============================================================================

subroutine ffdev_topology_EXP_update_nb_params(top)

    use ffdev_utils
    use ffdev_atomicdata

    implicit none
    type(TOPOLOGY)  :: top
    ! --------------------------------------------
    integer         :: i,ti,agti
    ! --------------------------------------------------------------------------

    ! first update PB from external sources
    select case(exp_pb_mode)
        case(EXP_PB_FREEOPT)
            ! nothing to do
    !---------------
        case(EXP_PB_DO)
            do i=1,top%nnb_types
                if( top%nb_types(i)%ti .ne. top%nb_types(i)%tj ) cycle
                ti   = top%nb_types(i)%ti
                agti = top%atom_types(ti)%glbtypeid
                top%nb_types(i)%pb = ffdev_atomicdata_do_bii(agti)
            end do
    !---------------
        case(EXP_PB_WO)
            do i=1,top%nnb_types
                if( top%nb_types(i)%ti .ne. top%nb_types(i)%tj ) cycle
                ti   = top%nb_types(i)%ti
                agti = top%atom_types(ti)%glbtypeid
                top%nb_types(i)%pb = ffdev_atomicdata_wo_bii(agti)
            end do
    !---------------
        case(EXP_PB_IP)
            do i=1,top%nnb_types
                if( top%nb_types(i)%ti .ne. top%nb_types(i)%tj ) cycle
                ti   = top%nb_types(i)%ti
                agti = top%atom_types(ti)%glbtypeid
                top%nb_types(i)%pb = ffdev_atomicdata_ip_bii(agti)
            end do
    !---------------
        case default
            call ffdev_utils_exit(DEV_ERR,1,'EXPPB mode not implemented in ffdev_topology_EXP_update_nb_params I!')
    end select

    ! apply combining rules if necessary
    select case(exp_pb_mode)
        case(EXP_PB_FREEOPT)
            if( ApplyCombiningRules ) then
                call ffdev_topology_EXP_update_nb_params_PB(top)
            end if
    !---------------
        case(EXP_PB_DO,EXP_PB_WO,EXP_PB_IP)
            if( .not. ApplyCombiningRules ) then
                ! we need to apply combining rules for unlike atoms
                call ffdev_utils_exit(DEV_ERR,1,'EXP_PB_DO, EXP_PB_WO, or EXP_PB_IP requires ApplyCombiningRules!')
            end if
            call ffdev_topology_EXP_update_nb_params_PB(top)
    !---------------
        case default
            call ffdev_utils_exit(DEV_ERR,1,'EXPPB mode not implemented in ffdev_topology_EXP_update_nb_params II!')
    end select

    ! and finally PA
    if( ApplyCombiningRules ) then
        call ffdev_topology_EXP_update_nb_params_PA(top)
    end if

end subroutine ffdev_topology_EXP_update_nb_params

!===============================================================================
! subroutine ffdev_topology_EXP_update_nb_params_PB
!===============================================================================

subroutine ffdev_topology_EXP_update_nb_params_PB(top)

    use ffdev_utils

    implicit none
    type(TOPOLOGY)  :: top
    ! --------------------------------------------
    integer         :: i,nbii,nbjj
    real(DEVDP)     :: pbii,pbjj,pbij
    ! --------------------------------------------------------------------------

   ! apply combining rules
    do i=1,top%nnb_types

        ! discard like atoms
        if( top%nb_types(i)%ti .eq. top%nb_types(i)%tj ) cycle

        ! get type parameters
        nbii = top%nb_types(i)%nbii
        nbjj = top%nb_types(i)%nbjj

        pbii  = top%nb_types(nbii)%pb
        pbjj  = top%nb_types(nbjj)%pb

        call ffdev_topology_EXP_apply_NB_comb_rules_PB(pbii,pbjj,pbij)

        top%nb_types(i)%pb = pbij

    end do

end subroutine ffdev_topology_EXP_update_nb_params_PB

!===============================================================================
! subroutine ffdev_topology_EXP_apply_NB_comb_rules_PB
!===============================================================================

subroutine ffdev_topology_EXP_apply_NB_comb_rules_PB(pbii,pbjj,pbij)

    use ffdev_utils

    implicit none
    real(DEVDP)     :: pbii,pbjj,pbij
    ! --------------------------------------------------------------------------

    select case(exp_comb_rules)
        case(EXP_COMB_RULE_AM)
            pbij = 0.5d0 * (pbii+pbjj)

        case(EXP_COMB_RULE_GS)
            if( pbii+pbjj .gt. 0 ) then
                pbij = 2.0d0 * pbii*pbjj/(pbii+pbjj)
            else
                pbij = 0.5d0 * (pbii+pbjj)
            end if

        case(EXP_COMB_RULE_BA)
            if( pbii+pbjj .gt. 0 ) then
                pbij = 2.0d0 * pbii*pbjj/(pbii+pbjj)
            else
                pbij = 0.5d0 * (pbii+pbjj)
            end if

        case(EXP_COMB_RULE_VS)
            pbij = sqrt(pbii*pbjj)

        case(EXP_COMB_RULE_D1,EXP_COMB_RULE_D2)
            if( (pbii**6 + pbjj**6) .gt. 0 ) then
                pbij = ( 2.0d0 * pbii**6 * pbjj**6 / (pbii**6 + pbjj**6) )**(1.0d0/6.0d0)
            else
                ! fallback
                pbij = 0.5d0 * (pbii+pbjj)
            end if

        case(EXP_COMB_RULE_W1,EXP_COMB_RULE_W2)
            if( (pbii**2 + pbjj**2) .gt. 0 ) then
                pbij = ( 2.0d0 * pbii**2 * pbjj**2 / (pbii**2 + pbjj**2) )**(1.0d0/2.0d0)
            else
                ! fallback
                pbij = 0.5d0 * (pbii+pbjj)
            end if

        case default
            call ffdev_utils_exit(DEV_ERR,1,'Not implemented nb_comb_rules in ffdev_topology_EXP_apply_NB_comb_rules_PB!')
    end select

end subroutine ffdev_topology_EXP_apply_NB_comb_rules_PB

!===============================================================================
! subroutine ffdev_topology_EXP_update_nb_params_PA
!===============================================================================

subroutine ffdev_topology_EXP_update_nb_params_PA(top)

    use ffdev_utils

    implicit none
    type(TOPOLOGY)  :: top
    ! --------------------------------------------
    integer         :: i,nbii,nbjj
    real(DEVDP)     :: pbii,pbjj,pbij,paii,pajj,paij,eaii,eajj,eaij
    ! --------------------------------------------------------------------------

   ! apply combining rules
    do i=1,top%nnb_types

        ! discard like atoms
        if( top%nb_types(i)%ti .eq. top%nb_types(i)%tj ) cycle

        ! get type parameters
        nbii = top%nb_types(i)%nbii
        nbjj = top%nb_types(i)%nbjj

        ! PA
        paii = top%nb_types(nbii)%pa
        pajj = top%nb_types(nbjj)%pa

        ! helpers for some rules
        pbii  = top%nb_types(nbii)%pb
        pbjj  = top%nb_types(nbjj)%pb
        pbij  = top%nb_types(i)%pb

        select case(exp_comb_rules)

            case(EXP_COMB_RULE_AM)
                paij = 0.5d0 * (paii+pajj)      ! paij is exponential of Aij, etc. ...

            case(EXP_COMB_RULE_GS)
                ! Z. Phys. D - Atoms, Molecules and Clusters 1, 91-101 (1986)
                eaii = exp(paii)
                eajj = exp(pajj)
                eaij = ( ( (eaii*pbii)**(1.0d0/pbii) * (eajj*pbjj)**(1.0d0/pbjj) )**(pbij/2.0d0) ) / pbij
                paij = log(eaij)                ! paij is exponential of Aij, etc. ...

            case(EXP_COMB_RULE_BA)
                ! Z. Phys. D - Atoms, Molecules and Clusters 1, 91-101 (1986)
                eaii = exp(paii)
                eajj = exp(pajj)
                eaij = ( (eaii)**(1.0d0/pbii) * (eajj)**(1.0d0/pbjj) )**(pbij/2.0d0)
                paij = log(eaij)                ! paij is exponential of Aij, etc. ...

            case(EXP_COMB_RULE_VS)
                ! DOI:10.1021/acs.jctc.6b00209J. Chem. Theory Comput.2016, 12, 3851−3870
                paij = 0.5d0 * (paii + pajj)    ! paij is exponential of Aij, etc. ...

            case(EXP_COMB_RULE_D1,EXP_COMB_RULE_W1)
                paij = 0.5d0 * (paii+pajj)      ! paij is exponential of Aij, etc. ...

            case(EXP_COMB_RULE_D2)
                paij = 0.5d0 * (paii+pajj) + 3.0d0/2.0d0 * ( log(pbii) + log(pbjj) ) - 3.0d0 * log(pbij)

            case(EXP_COMB_RULE_W2)
                paij = 0.5d0 * (paii+pajj) + 6.0d0/2.0d0 * ( log(pbii) + log(pbjj) ) - 6.0d0 * log(pbij)

            case default
                call ffdev_utils_exit(DEV_ERR,1,'Not implemented in ffdev_topology_EXP_apply_NB_comb_rules!')
        end select

        top%nb_types(i)%pa = paij

    end do

end subroutine ffdev_topology_EXP_update_nb_params_PA

! ------------------------------------------------------------------------------

end module ffdev_topology_exp

