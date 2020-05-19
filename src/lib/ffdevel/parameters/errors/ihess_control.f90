! ==============================================================================
! This file is part of FFDevel.
!    Copyright (C) 2019 Petr Kulhanek, kulhanek@chemi.muni.cz
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

module ffdev_err_ihess_control

use ffdev_constants
use ffdev_variables

contains

! ==============================================================================
! subroutine ffdev_err_ihess_ctrl
! ==============================================================================

subroutine ffdev_err_ihess_ctrl(fin)

    use ffdev_err_ihess_dat
    use ffdev_utils
    use prmfile

    implicit none
    type(PRMFILE_TYPE)  :: fin
    ! --------------------------------------------------------------------------

    write(DEV_OUT,*)
    write(DEV_OUT,10)

    if( .not. prmfile_open_section(fin,'ihess') ) then
        write(DEV_OUT,115) prmfile_onoff(EnableIHessError)
        write(DEV_OUT,135) prmfile_onoff(PrintIHessErrorSummary)
        write(DEV_OUT,125) IHessErrorsWeight
        write(DEV_OUT,145) prmfile_onoff(OnlyFFOptIHess)
        return
    end if

    if( prmfile_get_logical_by_key(fin,'enabled', EnableIHessError)) then
        write(DEV_OUT,110) prmfile_onoff(EnableIHessError)
    else
        write(DEV_OUT,115) prmfile_onoff(EnableIHessError)
    end if
    if( prmfile_get_logical_by_key(fin,'summary', PrintIHessErrorSummary)) then
        write(DEV_OUT,130) prmfile_onoff(PrintIHessErrorSummary)
    else
        write(DEV_OUT,135) prmfile_onoff(PrintIHessErrorSummary)
    end if
    if( prmfile_get_real8_by_key(fin,'weight', IHessErrorsWeight)) then
        write(DEV_OUT,120) IHessErrorsWeight
    else
        write(DEV_OUT,125) IHessErrorsWeight
    end if
    if( prmfile_get_logical_by_key(fin,'onlyffopt', OnlyFFOptIHess)) then
        write(DEV_OUT,140) prmfile_onoff(OnlyFFOptIHess)
    else
        write(DEV_OUT,145) prmfile_onoff(OnlyFFOptIHess)
    end if

 10 format('=== [ihess] ====================================================================')

110  format ('Bonds error (enabled)                  = ',a12)
115  format ('Bonds error (enabled)                  = ',a12,'                  (default)')
130  format ('Print ihess error summary (summary)    = ',a12)
135  format ('Print ihess error summary (summary)    = ',a12,'                  (default)')
120  format ('Bonds error weight (weight)            = ',f21.8)
125  format ('Bonds error weight (weight)            = ',f21.8,'         (default)')
140  format ('Only FFopt ihess (onlyffopt)           = ',a12)
145  format ('Only FFopt ihess (onlyffopt)           = ',a12,'                  (default)')

end subroutine ffdev_err_ihess_ctrl

! ------------------------------------------------------------------------------

end module ffdev_err_ihess_control
