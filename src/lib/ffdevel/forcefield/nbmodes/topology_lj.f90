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

module ffdev_topology_lj

use ffdev_topology_dat
use ffdev_constants
use ffdev_variables

contains

! ==============================================================================
! subroutine ffdev_topology_LJ_comb_rules_to_string
! ==============================================================================

character(80) function ffdev_topology_LJ_comb_rules_to_string(comb_rules)

    use ffdev_utils

    implicit none
    integer  :: comb_rules
    ! --------------------------------------------------------------------------

    select case(comb_rules)

        ! LJ potential
        case(LJ_COMB_RULE_LB)
            ffdev_topology_LJ_comb_rules_to_string = 'LB (Lorentz-Berthelot)'
        case(LJ_COMB_RULE_WH)
            ffdev_topology_LJ_comb_rules_to_string = 'WH (Waldman-Hagler)'
        case(LJ_COMB_RULE_KG)
            ffdev_topology_LJ_comb_rules_to_string = 'KG (Kong)'
        case(LJ_COMB_RULE_FH)
            ffdev_topology_LJ_comb_rules_to_string = 'FH (Fender-Halsey)'

        case default
            call ffdev_utils_exit(DEV_ERR,1,'Not implemented in ffdev_topology_LJ_comb_rules_to_string!')
    end select

end function ffdev_topology_LJ_comb_rules_to_string

! ==============================================================================
! function ffdev_topology_LJ_comb_rules_from_string
! ==============================================================================

integer function ffdev_topology_LJ_comb_rules_from_string(string)

    use ffdev_utils

    implicit none
    character(*)   :: string
    ! --------------------------------------------------------------------------

    select case(trim(string))
        case('LB')
            ffdev_topology_LJ_comb_rules_from_string = LJ_COMB_RULE_LB
        case('WH')
            ffdev_topology_LJ_comb_rules_from_string = LJ_COMB_RULE_WH
        case('KG')
            ffdev_topology_LJ_comb_rules_from_string = LJ_COMB_RULE_KG
        case('FH')
            ffdev_topology_LJ_comb_rules_from_string = LJ_COMB_RULE_FH

        case default
            call ffdev_utils_exit(DEV_ERR,1,'Not implemented "' // trim(string) //'" in ffdev_topology_LJ_comb_rules_from_string!')
    end select

end function ffdev_topology_LJ_comb_rules_from_string

! ==============================================================================
! subroutine ffdev_topology_LJ_alpha_mode_to_string
! ==============================================================================

character(80) function ffdev_topology_LJ_alpha_mode_to_string(alpha_mode)

    use ffdev_utils

    implicit none
    integer  :: alpha_mode
    ! --------------------------------------------------------------------------

    select case(alpha_mode)

        ! LJ potential
        case(LJ_ALPHA_CONST)
            ffdev_topology_LJ_alpha_mode_to_string = 'CONST'
        case(LJ_ALPHA_FREEOPT)
            ffdev_topology_LJ_alpha_mode_to_string = 'FREEOPT'
        case default
            call ffdev_utils_exit(DEV_ERR,1,'Not implemented in ffdev_topology_LJ_alpha_mode_to_string!')
    end select

end function ffdev_topology_LJ_alpha_mode_to_string

! ==============================================================================
! function ffdev_topology_LJ_alpha_mode_from_string
! ==============================================================================

integer function ffdev_topology_LJ_alpha_mode_from_string(string)

    use ffdev_utils

    implicit none
    character(*)   :: string
    ! --------------------------------------------------------------------------

    select case(trim(string))
        case('CONST')
            ffdev_topology_LJ_alpha_mode_from_string = LJ_ALPHA_CONST
        case('FREEOPT')
            ffdev_topology_LJ_alpha_mode_from_string = LJ_ALPHA_FREEOPT

        case default
            call ffdev_utils_exit(DEV_ERR,1,'Not implemented "' // trim(string) //'" in ffdev_topology_LJ_alpha_mode_from_string!')
    end select

end function ffdev_topology_LJ_alpha_mode_from_string

! ==============================================================================
! subroutine ffdev_topology_LJ_update_nb_params
! ==============================================================================

subroutine ffdev_topology_LJ_update_nb_params(top)

    use ffdev_utils

    implicit none
    type(TOPOLOGY)  :: top
    ! --------------------------------------------
    integer         :: i,nbii,nbjj
    real(DEVDP)     :: epsii,r0ii,epsjj,r0jj,epsij,r0ij,k,l
    real(DEVDP)     :: aii, ajj, aij
    ! --------------------------------------------------------------------------

    if( .not. ApplyCombiningRules ) return

    ! apply combining rules
    do i=1,top%nnb_types

        ! discard like atoms
        if( top%nb_types(i)%ti .eq. top%nb_types(i)%tj ) cycle

        ! get type parameters
        nbii = top%nb_types(i)%nbii
        nbjj = top%nb_types(i)%nbjj

        r0ii  = top%nb_types(nbii)%r0
        epsii = top%nb_types(nbii)%eps
        aii   = top%nb_types(nbii)%alpha

        r0jj  = top%nb_types(nbjj)%r0
        epsjj = top%nb_types(nbjj)%eps
        ajj   = top%nb_types(nbjj)%alpha

        select case(lj_comb_rules)
            case(LJ_COMB_RULE_LB)
                r0ij = (r0ii+r0jj)*0.5d0
                epsij = sqrt(epsii*epsjj)
            case(LJ_COMB_RULE_WH)
                r0ij = ((r0ii**6 + r0jj**6)*0.5d0)**(1.0d0/6.0d0)
                epsij = sqrt( epsii*r0ii**6 * epsjj*r0jj**6 )/r0ij**6
            case(LJ_COMB_RULE_KG)
                k = sqrt(epsii*r0ii**6 * epsjj*r0jj**6)
                l = ( ( (epsii*r0ii**12)**(1.0d0/13.0d0) + (epsjj*r0jj**12)**(1.0d0/13.0d0) )*0.5d0 )**13
                r0ij = (l/k)**(1.0d0/6.0d0)
                epsij = k / (r0ij**6)
            case(LJ_COMB_RULE_FH)
                r0ij = (r0ii+r0jj)*0.5d0
                epsij = 2.0d0*epsii*epsjj/(epsii+epsjj)
            case default
                call ffdev_utils_exit(DEV_ERR,1,'Not implemented in ffdev_topology_LJ_update_nb_params!')
        end select

        ! only for probing in exp-6 mode
        aij = 0.5d0 * (aii + ajj)

        top%nb_types(i)%r0    = r0ij
        top%nb_types(i)%eps   = epsij
        top%nb_types(i)%alpha = aij
    end do

end subroutine ffdev_topology_LJ_update_nb_params

! ------------------------------------------------------------------------------

end module ffdev_topology_lj

