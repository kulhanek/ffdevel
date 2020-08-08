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

module ffdev_topology_pen

use ffdev_topology_dat
use ffdev_constants
use ffdev_variables

contains

! ==============================================================================
! subroutine ffdev_topology_PEN_pa_mode_to_string
! ==============================================================================

character(80) function ffdev_topology_PEN_pa_mode_to_string(lpen_mode)

    use ffdev_utils

    implicit none
    integer  :: lpen_mode
    ! --------------------------------------------------------------------------

    select case(lpen_mode)
        case(PEN_PA_FREEOPT)
            ffdev_topology_PEN_pa_mode_to_string = 'FREEOPT - use optimized PA per type'
        case(PEN_PA_CONST)
            ffdev_topology_PEN_pa_mode_to_string = 'CONST - set to pen_fa'
        case(PEN_PA_ADBII)
            ffdev_topology_PEN_pa_mode_to_string = 'ADBII - Atomic database Bii'
        case(PEN_PA_ADRCII)
            ffdev_topology_PEN_pa_mode_to_string = 'ADRCII - Atomic database Rcii'
        case(PEN_PA_COUPLED)
            ffdev_topology_PEN_pa_mode_to_string = 'COUPLED - PA coupled by pen_fa to EXP(PB)'
        case default
            call ffdev_utils_exit(DEV_ERR,1,'Not implemented in ffdev_topology_PEN_pa_mode_to_string!')
    end select

end function ffdev_topology_PEN_pa_mode_to_string

! ==============================================================================
! subroutine ffdev_topology_PEN_pa_mode_from_string
! ==============================================================================

integer function ffdev_topology_PEN_pa_mode_from_string(string)

    use ffdev_utils

    implicit none
    character(*)   :: string
    ! --------------------------------------------------------------------------

    select case(trim(string))
        case('FREEOPT')
            ffdev_topology_PEN_pa_mode_from_string = PEN_PA_FREEOPT
        case('CONST')
            ffdev_topology_PEN_pa_mode_from_string = PEN_PA_CONST
        case('ADBII')
            ffdev_topology_PEN_pa_mode_from_string = PEN_PA_ADBII
        case('ADRCII')
            ffdev_topology_PEN_pa_mode_from_string = PEN_PA_ADRCII
        case('COUPLED')
            ffdev_topology_PEN_pa_mode_from_string = PEN_PA_COUPLED
        case default
            call ffdev_utils_exit(DEV_ERR,1,'Not implemented "' // trim(string) //'" in ffdev_topology_PEN_pa_mode_from_string!')
    end select

end function ffdev_topology_PEN_pa_mode_from_string

! ==============================================================================
! subroutine ffdev_topology_PEN_pb_mode_to_string
! ==============================================================================

character(80) function ffdev_topology_PEN_pb_mode_to_string(lpen_mode)

    use ffdev_utils

    implicit none
    integer  :: lpen_mode
    ! --------------------------------------------------------------------------

    select case(lpen_mode)
        case(PEN_PB_FREEOPT)
            ffdev_topology_PEN_pb_mode_to_string = 'FREEOPT - use optimized PB per type'
        case(PEN_PB_CONST)
            ffdev_topology_PEN_pb_mode_to_string = 'CONST - set to pen_fb'
        case(PEN_PB_COUPLED)
            ffdev_topology_PEN_pb_mode_to_string = 'COUPLED - PB coupled by pen_fb to PA'
        case default
            call ffdev_utils_exit(DEV_ERR,1,'Not implemented in ffdev_topology_PEN_pb_mode_to_string!')
    end select

end function ffdev_topology_PEN_pb_mode_to_string

! ==============================================================================
! subroutine ffdev_topology_PEN_pb_mode_from_string
! ==============================================================================

integer function ffdev_topology_PEN_pb_mode_from_string(string)

    use ffdev_utils

    implicit none
    character(*)   :: string
    ! --------------------------------------------------------------------------

    select case(trim(string))
        case('FREEOPT')
            ffdev_topology_PEN_pb_mode_from_string = PEN_PB_FREEOPT
        case('CONST')
            ffdev_topology_PEN_pb_mode_from_string = PEN_PB_CONST
        case('COUPLED')
            ffdev_topology_PEN_pb_mode_from_string = PEN_PB_COUPLED
        case default
            call ffdev_utils_exit(DEV_ERR,1,'Not implemented "' // trim(string) //'" in ffdev_topology_PEN_pb_mode_from_string!')
    end select

end function ffdev_topology_PEN_pb_mode_from_string

! ------------------------------------------------------------------------------

integer function ffdev_topology_PEN_get_valZ(gti)

    use ffdev_utils
    use ffdev_parameters_dat

    implicit none
    integer     :: gti
    ! --------------------------------------------
    integer     :: z
    ! --------------------------------------------------------------------------

    ffdev_topology_PEN_get_valZ = 1 ! default to H

    ! get Z
    z = types(gti)%z
    if( z .le. 0 ) then
        call ffdev_utils_exit(DEV_ERR,1,'Z is out-of-range in ffdev_topology_PEN_get_valZ')
    end if
    if( z .le. 2 ) then
        ! H-He
        ffdev_topology_PEN_get_valZ = z
        return
    end if
    if( z .le. 10 ) then
        ! Li-Ne
        ffdev_topology_PEN_get_valZ = z - 2
        return
    end if
    if( z .le. 18 ) then
        ! Na-Ar
        ffdev_topology_PEN_get_valZ = z - 10
        return
    end if

    call ffdev_utils_exit(DEV_ERR,1,'Not implemented in ffdev_topology_PEN_get_valZ')

end function ffdev_topology_PEN_get_valZ

! ==============================================================================
! subroutine ffdev_topology_PEN_update_nb_params
! ==============================================================================

subroutine ffdev_topology_PEN_update_nb_params(top)

    use ffdev_utils
    use ffdev_atomicdata
    use ffdev_disp_dat
    use ffdev_topology_utils

    implicit none
    type(TOPOLOGY)  :: top
    ! --------------------------------------------
    integer         :: i,nbtii
    real(DEVDP)     :: r
    ! --------------------------------------------------------------------------

    do i=1,top%natom_types
        select case(pen_pa_mode)
            case(PEN_PA_FREEOPT)
                ! nothing to do
            case(PEN_PA_CONST)
                top%atom_types(i)%pen_pa = pen_fa
            case(PEN_PA_ADBII)
                top%atom_types(i)%pen_pa = pen_fa * ffdev_atomicdata_bii(top%atom_types(i)%glbtypeid)
            case(PEN_PA_ADRCII)
                r = ffdev_atomicdata_rcii(top%atom_types(i)%glbtypeid,damp_fa,damp_fb)
                if( r .le. 0.0d0 ) r = 1.0d0
                top%atom_types(i)%pen_pa = pen_fa / r
            case(PEN_PA_COUPLED)
                nbtii = ffdev_topology_find_nbtype_by_tindex(top,i,i)
                top%atom_types(i)%pen_pa = pen_fa * top%nb_types(nbtii)%PB
            case default
                call ffdev_utils_exit(DEV_ERR,1,'PEN_PA_MODE mode not implemented in ffdev_topology_PEN_update_nb_params!')
        end select
    ! --------
        select case(pen_pb_mode)
            case(PEN_PB_FREEOPT)
                ! nothing to do
            case(PEN_PB_CONST)
                top%atom_types(i)%pen_pb = pen_fb
            case(PEN_PB_COUPLED)
                top%atom_types(i)%pen_pb = pen_fb * top%atom_types(i)%pen_pa
            case default
                call ffdev_utils_exit(DEV_ERR,1,'PEN_PB_MODE mode not implemented in ffdev_topology_PEN_update_nb_params!')
        end select
    end do

end subroutine ffdev_topology_PEN_update_nb_params

! ------------------------------------------------------------------------------

end module ffdev_topology_pen

