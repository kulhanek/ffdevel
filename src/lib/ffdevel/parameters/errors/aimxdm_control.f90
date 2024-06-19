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

module ffdev_err_aimxdm_control

use ffdev_constants
use ffdev_variables

contains

! ==============================================================================
! subroutine ffdev_err_aimxdm_ctrl
! ==============================================================================

subroutine ffdev_err_aimxdm_ctrl(fin)

    use ffdev_err_aimxdm_dat
    use ffdev_errors_dat
    use ffdev_utils
    use ffdev_errors_utils
    use prmfile

    implicit none
    type(PRMFILE_TYPE)          :: fin
    ! --------------------------------------------------------------------------

    write(DEV_OUT,*)
    write(DEV_OUT,10)

    if( .not. prmfile_open_section(fin,'aimxdm') ) then
        write(DEV_OUT,115) prmfile_onoff(EnableAIMXDMError)
        write(DEV_OUT,135) prmfile_onoff(PrintAIMXDMErrorSummary)
        write(DEV_OUT,125) AIMXDMErrorWeight
        write(DEV_OUT,145) prmfile_onoff(AIMXDMBurriedOnly)
        return
    end if

    if( prmfile_get_logical_by_key(fin,'enabled', EnableAIMXDMError)) then
        write(DEV_OUT,110) prmfile_onoff(EnableAIMXDMError)
    else
        write(DEV_OUT,115) prmfile_onoff(EnableAIMXDMError)
    end if

    if( prmfile_get_logical_by_key(fin,'summary', PrintAIMXDMErrorSummary)) then
        write(DEV_OUT,130) prmfile_onoff(PrintAIMXDMErrorSummary)
    else
        write(DEV_OUT,135) prmfile_onoff(PrintAIMXDMErrorSummary)
    end if

    if( prmfile_get_real8_by_key(fin,'weight', AIMXDMErrorWeight)) then
        write(DEV_OUT,120) AIMXDMErrorWeight
    else
        write(DEV_OUT,125) AIMXDMErrorWeight
    end if

    if( prmfile_get_logical_by_key(fin,'buried', AIMXDMBurriedOnly)) then
        write(DEV_OUT,140) prmfile_onoff(AIMXDMBurriedOnly)
    else
        write(DEV_OUT,145) prmfile_onoff(AIMXDMBurriedOnly)
    end if

 10 format('=== [aimxdm] ===================================================================')

110  format ('LJ R0 penalty (enabled)                = ',a12)
115  format ('LJ R0 penalty (enabled)                = ',a12,'                  (default)')
130  format ('LJ R0 penalty summary (summary)        = ',a12)
135  format ('LJ R0 penalty summary (summary)        = ',a12,'                  (default)')

120  format ('LJ R0 penalty weight (weight)          = ',f21.8)
125  format ('LJ R0 penalty weight (weight)          = ',f21.8,'         (default)')

140  format ('LJ R0 penalty only buried (buried)     = ',a12)
145  format ('LJ R0 penalty only buried (buried)     = ',a12,'                  (default)')

end subroutine ffdev_err_aimxdm_ctrl

! ------------------------------------------------------------------------------

end module ffdev_err_aimxdm_control
