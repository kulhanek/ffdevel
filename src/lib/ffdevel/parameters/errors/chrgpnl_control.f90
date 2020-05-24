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

module ffdev_err_chrgpnl_control

use ffdev_constants
use ffdev_variables

contains

! ==============================================================================
! subroutine ffdev_err_chrgpnl_ctrl
! ==============================================================================

subroutine ffdev_err_chrgpnl_ctrl(fin)

    use ffdev_err_chrgpnl_dat
    use ffdev_errors_dat
    use ffdev_utils
    use ffdev_errors_utils
    use prmfile

    implicit none
    type(PRMFILE_TYPE)          :: fin
    ! --------------------------------------------------------------------------

    write(DEV_OUT,*)
    write(DEV_OUT,10)

    if( .not. prmfile_open_section(fin,'chrgpnl') ) then
        write(DEV_OUT,115) prmfile_onoff(EnableChrgPnlError)
        write(DEV_OUT,135) prmfile_onoff(PrintChrgPnlErrorSummary)

        write(DEV_OUT,125) ChrgPnlErrorWeight
        return
    end if

    if( prmfile_get_logical_by_key(fin,'enabled', EnableChrgPnlError)) then
        write(DEV_OUT,110) prmfile_onoff(EnableChrgPnlError)
    else
        write(DEV_OUT,115) prmfile_onoff(EnableChrgPnlError)
    end if
    if( prmfile_get_logical_by_key(fin,'summary', PrintChrgPnlErrorSummary)) then
        write(DEV_OUT,130) prmfile_onoff(PrintChrgPnlErrorSummary)
    else
        write(DEV_OUT,135) prmfile_onoff(PrintChrgPnlErrorSummary)
    end if

    if( prmfile_get_real8_by_key(fin,'weight', ChrgPnlErrorWeight)) then
        write(DEV_OUT,120) ChrgPnlErrorWeight
    else
        write(DEV_OUT,125) ChrgPnlErrorWeight
    end if

 10 format('=== [chrgpnl] =================================================================')

110  format ('Charge penalty (enabled)               = ',a12)
115  format ('Charge penalty (enabled)               = ',a12,'                  (default)')
130  format ('Charge penalty summary (summary)       = ',a12)
135  format ('Charge penalty summary (summary)       = ',a12,'                  (default)')

120  format ('Charge penalty weight (weight)         = ',f21.8)
125  format ('Charge penalty weight (weight)         = ',f21.8,'         (default)')

end subroutine ffdev_err_chrgpnl_ctrl

! ------------------------------------------------------------------------------

end module ffdev_err_chrgpnl_control
