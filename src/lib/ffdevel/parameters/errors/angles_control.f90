! ==============================================================================
! This file is part of FFDevel.
!    Copyright (C) 2018 Petr Kulhanek, kulhanek@chemi.muni.cz
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

module ffdev_err_angles_control

use ffdev_constants

contains

! ==============================================================================
! subroutine ffdev_err_angles_ctrl
! ==============================================================================

subroutine ffdev_err_angles_ctrl(fin)

    use ffdev_err_angles_dat
    use ffdev_utils
    use prmfile

    implicit none
    type(PRMFILE_TYPE)  :: fin
    ! --------------------------------------------------------------------------

    write(DEV_OUT,*)
    write(DEV_OUT,10)

    if( .not. prmfile_open_section(fin,'angles') ) then
        write(DEV_OUT,115) prmfile_onoff(EnableAngleError)
        write(DEV_OUT,135) prmfile_onoff(PrintAngleErrorSummary)
        write(DEV_OUT,125) AngleErrorWeight
        return
    end if

    if( prmfile_get_logical_by_key(fin,'enabled', EnableAngleError)) then
        write(DEV_OUT,110) prmfile_onoff(EnableAngleError)
    else
        write(DEV_OUT,115) prmfile_onoff(EnableAngleError)
    end if
    if( prmfile_get_logical_by_key(fin,'summary', PrintAngleErrorSummary)) then
        write(DEV_OUT,130) prmfile_onoff(PrintAngleErrorSummary)
    else
        write(DEV_OUT,135) prmfile_onoff(PrintAngleErrorSummary)
    end if
    if( prmfile_get_real8_by_key(fin,'weight', AngleErrorWeight)) then
        write(DEV_OUT,120) AngleErrorWeight
    else
        write(DEV_OUT,125) AngleErrorWeight
    end if

 10 format('=== [angles] ===================================================================')

110  format ('Angle error (enabled)                  = ',a12)
115  format ('Angle error (enabled)                  = ',a12,'                  (default)')
130  format ('Print angle error summary (summary)    = ',a12)
135  format ('Print angle error summary (summary)    = ',a12,'                  (default)')
120  format ('Angle error weight (weight)            = ',f21.8)
125  format ('Angle error weight (weight)            = ',f21.8,'         (default)')

end subroutine ffdev_err_angles_ctrl

! ------------------------------------------------------------------------------

end module ffdev_err_angles_control
