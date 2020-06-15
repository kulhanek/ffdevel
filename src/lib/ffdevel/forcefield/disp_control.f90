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

module ffdev_disp_control

use ffdev_sizes
use ffdev_constants
use ffdev_variables

contains

! ==============================================================================
! subroutine ffdev_disp_ctrl
! ==============================================================================

subroutine ffdev_disp_ctrl(fin,exec)

    use prmfile
    use ffdev_utils
    use ffdev_xdm_dat
    use ffdev_mmd3_dat
    use ffdev_disp_dat
    use ffdev_disp
    use ffdev_parameters_dat

    implicit none
    type(PRMFILE_TYPE)          :: fin
    logical                     :: exec
    ! --------------------------------------------
    integer                     :: alloc_stat
    character(PRMFILE_MAX_PATH) :: string
    ! --------------------------------------------------------------------------

    write(DEV_OUT,*)
    write(DEV_OUT,10)

    if( .not. prmfile_open_section(fin,'disp') ) then
        write(DEV_OUT,25) ffdev_disp_cxsource_to_string(cx_source)
        write(DEV_OUT,35) ffdev_disp_rcsource_to_string(rc_source)
        return
    end if

    if( prmfile_get_string_by_key(fin,'cx_source', string)) then
        cx_source = ffdev_disp_cxsource_from_string(string)
        write(DEV_OUT,20) ffdev_disp_cxsource_to_string(cx_source)
    else
        write(DEV_OUT,25) ffdev_disp_cxsource_to_string(cx_source)
    end if

    if( prmfile_get_string_by_key(fin,'rc_source', string)) then
        rc_source = ffdev_disp_rcsource_from_string(string)
        write(DEV_OUT,30) ffdev_disp_rcsource_to_string(rc_source)
    else
        write(DEV_OUT,35) ffdev_disp_rcsource_to_string(rc_source)
    end if

    if( .not. exec ) return

! execute
    if( .not. disp_data_loaded ) then
        allocate( disp_pairs(ntypes,ntypes), stat = alloc_stat)
        if(alloc_stat .ne. 0) then
            call ffdev_utils_exit(DEV_ERR,1,'Unable to allocate memory for disp_pairs in ffdev_parameters_ctrl_disp!')
        end if
        disp_data_loaded = .true.
    end if

    disp_pairs(:,:)%c6 = 0.0d0
    disp_pairs(:,:)%c8 = 0.0d0
    disp_pairs(:,:)%c10 = 0.0d0
    disp_pairs(:,:)%rc = 0.0d0

    select case(cx_source)
        case(NB_CX_XDM)
            disp_pairs(:,:)%c6  = xdm_pairs(:,:)%c6ave
            disp_pairs(:,:)%c8  = xdm_pairs(:,:)%c8ave
            disp_pairs(:,:)%c10 = xdm_pairs(:,:)%c10ave
        case(NB_CX_MMD3)
            disp_pairs(:,:)%c6  = mmd3_pairs(:,:)%c6ave
            disp_pairs(:,:)%c8  = mmd3_pairs(:,:)%c8ave
        case default
            call ffdev_utils_exit(DEV_ERR,1,'Not implemented in ffdev_parameters_ctrl_disp I!')
    end select

    select case(rc_source)
        case(NB_RC_XDM)
            disp_pairs(:,:)%Rc  = xdm_pairs(:,:)%Rc
        case(NB_RC_XDM_POL)
            disp_pairs(:,:)%Rc  = xdm_pairs(:,:)%Rvdw
        case(NB_RC_XDM_VOL)
            disp_pairs(:,:)%Rc  = xdm_pairs(:,:)%Rvol
        case(NB_RC_MMD3)
            disp_pairs(:,:)%Rc  = mmd3_pairs(:,:)%Rc
        case default
            call ffdev_utils_exit(DEV_ERR,1,'Not implemented in ffdev_parameters_ctrl_disp II!')
    end select

    return

 10 format('=== [disp] =====================================================================')

 20  format ('Cx source (cx_source)                    = ',a12)
 25  format ('Cx source (cx_source)                    = ',a12,'                (default)')

 30  format ('Rc source (rc_source)                    = ',a12)
 35  format ('Rc source (rc_source)                    = ',a12,'                (default)')

end subroutine ffdev_disp_ctrl

! ==============================================================================

end module ffdev_disp_control

