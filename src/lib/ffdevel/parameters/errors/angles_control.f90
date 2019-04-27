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
        write(DEV_OUT,115) prmfile_onoff(EnableAnglesError)
        write(DEV_OUT,135) prmfile_onoff(PrintAnglesErrorSummary)
        write(DEV_OUT,125) AngleErrorsWeight
        write(DEV_OUT,145) prmfile_onoff(OnlyFFOptAngles)
        return
    end if

    if( prmfile_get_logical_by_key(fin,'enabled', EnableAnglesError)) then
        write(DEV_OUT,110) prmfile_onoff(EnableAnglesError)
    else
        write(DEV_OUT,115) prmfile_onoff(EnableAnglesError)
    end if
    if( prmfile_get_logical_by_key(fin,'summary', PrintAnglesErrorSummary)) then
        write(DEV_OUT,130) prmfile_onoff(PrintAnglesErrorSummary)
    else
        write(DEV_OUT,135) prmfile_onoff(PrintAnglesErrorSummary)
    end if
    if( prmfile_get_real8_by_key(fin,'weight', AngleErrorsWeight)) then
        write(DEV_OUT,120) AngleErrorsWeight
    else
        write(DEV_OUT,125) AngleErrorsWeight
    end if
    if( prmfile_get_logical_by_key(fin,'onlyffopt', OnlyFFOptAngles)) then
        write(DEV_OUT,140) prmfile_onoff(OnlyFFOptAngles)
    else
        write(DEV_OUT,145) prmfile_onoff(OnlyFFOptAngles)
    end if

 10 format('=== [angles] ===================================================================')

110  format ('Angle error (enabled)                  = ',a12)
115  format ('Angle error (enabled)                  = ',a12,'                  (default)')
130  format ('Print angle error summary (summary)    = ',a12)
135  format ('Print angle error summary (summary)    = ',a12,'                  (default)')
120  format ('Angle error weight (weight)            = ',f21.8)
125  format ('Angle error weight (weight)            = ',f21.8,'         (default)')
140  format ('Only FFopt angles (onlyffopt)          = ',a12)
145  format ('Only FFopt angles (onlyffopt)          = ',a12,'                  (default)')

end subroutine ffdev_err_angles_ctrl

! ------------------------------------------------------------------------------

end module ffdev_err_angles_control
