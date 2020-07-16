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
        case(EXP_PB_DO_FULL)
            ffdev_topology_EXP_pb_mode_to_string = 'DO-FULL - PB from density overlap including unlike'
        case(EXP_PB_IP)
            ffdev_topology_EXP_pb_mode_to_string = 'IP - Derived from ionization potentials'
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
        case('DO-FULL')
            ffdev_topology_EXP_pb_mode_from_string = EXP_PB_DO_FULL
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
        case(REP_COMB_RULE_BDK)
            ffdev_topology_EXP_comb_rules_to_string = 'BDK (Bouchal-Durnik-Kulhanek)'
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
        case('BDK')
            ffdev_topology_EXP_comb_rules_from_string = REP_COMB_RULE_BDK
        case default
            call ffdev_utils_exit(DEV_ERR,1,'Not implemented "' // trim(string) //'" in ffdev_topology_EXP_comb_rules_from_string!')
    end select

end function ffdev_topology_EXP_comb_rules_from_string

!! ==============================================================================
!! subroutine ffdev_topology_apply_NB_comb_rules_EXP
!! ==============================================================================
!
!subroutine ffdev_topology_apply_NB_comb_rules_EXPTT(top,comb_rules)
!
!    use ffdev_utils
!
!    implicit none
!    type(TOPOLOGY)  :: top
!    integer         :: comb_rules
!    ! --------------------------------------------
!    integer         :: i,nbii,nbjj
!    real(DEVDP)     :: paii,paij,pajj
!    real(DEVDP)     :: pbii,pbij,pbjj
!    real(DEVDP)     :: tbii,tbij,tbjj
!    real(DEVDP)     :: eaii,eaij,eajj
!    ! --------------------------------------------------------------------------
!
!    ! apply combining rules
!    do i=1,top%nnb_types
!        if( top%nb_types(i)%ti .ne. top%nb_types(i)%tj ) then
!
!            ! get type parameters
!            nbii = ffdev_topology_find_nbtype_by_tindex(top,top%nb_types(i)%ti,top%nb_types(i)%ti)
!            nbjj = ffdev_topology_find_nbtype_by_tindex(top,top%nb_types(i)%tj,top%nb_types(i)%tj)
!
!            paii = top%nb_types(nbii)%pa
!            pbii = top%nb_types(nbii)%pb
!            tbii = top%nb_types(nbii)%tb
!
!            pajj = top%nb_types(nbjj)%pa
!            pbjj = top%nb_types(nbjj)%pb
!            tbjj = top%nb_types(nbjj)%tb
!
!            select case(comb_rules)
!
!                case(COMB_RULE_AM)
!                    paij = 0.5d0 * (paii+pajj)
!                    pbij = 0.5d0 * (pbii+pbjj)
!                    tbij = 0.5d0 * (tbii+tbjj)
!
!                case(COMB_RULE_GS)
!                    if( pbii+pbjj .gt. 0 ) then
!                        pbij = 2.0d0 * pbii*pbjj/(pbii+pbjj)
!                    else
!                        pbij = 0.5d0 * (pbii+pbjj)
!                    end if
!
!                    ! Z. Phys. D - Atoms, Molecules and Clusters 1, 91-101 (1986)
!                    eaii = exp(paii)
!                    eajj = exp(pajj)
!                    eaij = ( ( (eaii*pbii)**(1.0d0/pbii) * (eajj*pbjj)**(1.0d0/pbjj) )**(pbij/2.0d0) ) / pbij
!                    paij = log(eaij)
!
!                    if( tbii+tbjj .gt. 0 ) then
!                        tbij = 2.0d0 * tbii*tbjj/(tbii+tbjj)
!                    else
!                        tbij = 0.5d0 * (tbii+tbjj)
!                    end if
!
!                case(COMB_RULE_BA)
!                    if( pbii+pbjj .gt. 0 ) then
!                        pbij = 2.0d0 * pbii*pbjj/(pbii+pbjj)
!                    else
!                        pbij = 0.5d0 * (pbii+pbjj)
!                    end if
!
!                    ! Z. Phys. D - Atoms, Molecules and Clusters 1, 91-101 (1986)
!                    eaii = exp(paii)
!                    eajj = exp(pajj)
!                    eaij = ( (eaii)**(1.0d0/pbii) * (eajj)**(1.0d0/pbjj) )**(pbij/2.0d0)
!                    paij = log(eaij)
!
!                    if( tbii+tbjj .gt. 0 ) then
!                        tbij = 2.0d0 * tbii*tbjj/(tbii+tbjj)
!                    else
!                        tbij = 0.5d0 * (tbii+tbjj)
!                    end if
!
!                case(COMB_RULE_VS)
!                    ! DOI:10.1021/acs.jctc.6b00209J. Chem. Theory Comput.2016, 12, 3851−3870
!                    paij = 0.5d0 * (paii + pajj)  ! paij is exponential of Aij, etc. ...
!                    pbij = sqrt(pbii*pbjj)
!                    tbij = sqrt(tbii*tbjj)
!
!                case default
!                    call ffdev_utils_exit(DEV_ERR,1,'Not implemented in ffdev_topology_apply_NB_comb_rules_EXPTT!')
!            end select
!
!            top%nb_types(i)%pa = paij
!            top%nb_types(i)%pb = pbij
!            top%nb_types(i)%tb = tbij
!
!        end if
!    end do
!
!end subroutine ffdev_topology_apply_NB_comb_rules_EXPTT

!===============================================================================
! subroutine ffdev_topology_EXP_apply_NB_comb_rules_PB
!===============================================================================

subroutine ffdev_topology_EXP_apply_NB_comb_rules_PB(pbii,pbjj,pbij)

    use ffdev_utils

    implicit none
    real(DEVDP)     :: pbii,pbjj,pbij
    ! --------------------------------------------------------------------------

    !write(*,*) pbii,pbjj
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

        case default
            call ffdev_utils_exit(DEV_ERR,1,'Not implemented nb_comb_rules in ffdev_topology_EXP_apply_NB_comb_rules_PB!')
    end select

end subroutine ffdev_topology_EXP_apply_NB_comb_rules_PB

! ------------------------------------------------------------------------------

end module ffdev_topology_exp

