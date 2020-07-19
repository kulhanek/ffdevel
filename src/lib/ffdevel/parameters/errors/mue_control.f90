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

module ffdev_err_mue_control

use ffdev_constants
use ffdev_variables

contains

! ==============================================================================
! subroutine ffdev_err_mue_ctrl
! ==============================================================================

subroutine ffdev_err_mue_ctrl(fin)

    use ffdev_err_mue_dat
    use ffdev_errors_dat
    use ffdev_utils
    use ffdev_errors_utils
    use prmfile

    implicit none
    type(PRMFILE_TYPE)          :: fin
    ! --------------------------------------------
    character(PRMFILE_MAX_PATH) :: string
    ! --------------------------------------------------------------------------

    write(DEV_OUT,*)
    write(DEV_OUT,10)

    if( .not. prmfile_open_section(fin,'mue') ) then
        write(DEV_OUT,115) prmfile_onoff(EnableMUEError)
        write(DEV_OUT,125) MUEErrorWeight
    end if

    if( prmfile_get_logical_by_key(fin,'enabled', EnableMUEError)) then
        write(DEV_OUT,110) prmfile_onoff(EnableMUEError)
    else
        write(DEV_OUT,115) prmfile_onoff(EnableMUEError)
    end if
    if( prmfile_get_real8_by_key(fin,'weight', MUEErrorWeight)) then
        write(DEV_OUT,120) MUEErrorWeight
    else
        write(DEV_OUT,125) MUEErrorWeight
    end if

 10 format('=== [mue] =================================================================')

110  format ('Energy error (enabled)                 = ',a12)
115  format ('Energy error (enabled)                 = ',a12,'                  (default)')
120  format ('Energy error weight (weight)           = ',f21.8)
125  format ('Energy error weight (weight)           = ',f21.8,'         (default)')

end subroutine ffdev_err_mue_ctrl

! ------------------------------------------------------------------------------

end module ffdev_err_mue_control
