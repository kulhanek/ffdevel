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
    use ffdev_disp_dat
    use ffdev_disp

    implicit none
    type(PRMFILE_TYPE)          :: fin
    logical                     :: exec
    ! --------------------------------------------

    character(PRMFILE_MAX_PATH) :: string
    ! --------------------------------------------------------------------------

    write(DEV_OUT,*)
    write(DEV_OUT,10)

    if( .not. prmfile_open_section(fin,'disp') ) then
        write(DEV_OUT,25) ffdev_disp_cxsource_to_string(cx_source)
        write(DEV_OUT,35) ffdev_disp_rcsource_to_string(rc_source)
        call ffdev_disp_update_db
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

    call ffdev_disp_update_db

    return

 10 format('=== [disp] =====================================================================')

 20  format ('Cx source (cx_source)                    = ',a12)
 25  format ('Cx source (cx_source)                    = ',a12,'                (default)')

 30  format ('Rc source (rc_source)                    = ',a12)
 35  format ('Rc source (rc_source)                    = ',a12,'                (default)')

end subroutine ffdev_disp_ctrl

! ==============================================================================

end module ffdev_disp_control

