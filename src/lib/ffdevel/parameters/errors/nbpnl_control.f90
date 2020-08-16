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

module ffdev_err_nbpnl_control

use ffdev_constants
use ffdev_variables

contains

! ==============================================================================
! subroutine ffdev_err_nbpnl_ctrl
! ==============================================================================

subroutine ffdev_err_nbpnl_ctrl(fin)

    use ffdev_err_nbpnl_dat
    use ffdev_errors_dat
    use ffdev_utils
    use ffdev_errors_utils
    use prmfile

    implicit none
    type(PRMFILE_TYPE)          :: fin
    ! --------------------------------------------------------------------------

    write(DEV_OUT,*)
    write(DEV_OUT,10)

    if( .not. prmfile_open_section(fin,'nbpnl') ) then
        write(DEV_OUT,115) prmfile_onoff(EnableNBPnlError)
        write(DEV_OUT,135) prmfile_onoff(PrintNBPnlErrorSummary)
        write(DEV_OUT,145) NBPnlErrorTempFactor
        write(DEV_OUT,125) NBPnlErrorWeight
        return
    end if

    if( prmfile_get_logical_by_key(fin,'enabled', EnableNBPnlError)) then
        write(DEV_OUT,110) prmfile_onoff(EnableNBPnlError)
    else
        write(DEV_OUT,115) prmfile_onoff(EnableNBPnlError)
    end if

    if( prmfile_get_logical_by_key(fin,'summary', PrintNBPnlErrorSummary)) then
        write(DEV_OUT,130) prmfile_onoff(PrintNBPnlErrorSummary)
    else
        write(DEV_OUT,135) prmfile_onoff(PrintNBPnlErrorSummary)
    end if

    if( prmfile_get_real8_by_key(fin,'temp', NBPnlErrorTempFactor)) then
        write(DEV_OUT,140) NBPnlErrorTempFactor
    else
        write(DEV_OUT,145) NBPnlErrorTempFactor
    end if

    if( prmfile_get_real8_by_key(fin,'weight', NBPnlErrorWeight)) then
        write(DEV_OUT,120) NBPnlErrorWeight
    else
        write(DEV_OUT,125) NBPnlErrorWeight
    end if

 10 format('=== [nbpnl] ====================================================================')

110  format ('QNB penalty (enabled)                  = ',a12)
115  format ('QNB penalty (enabled)                  = ',a12,'                  (default)')
130  format ('QNB penalty summary (summary)          = ',a12)
135  format ('QNB penalty summary (summary)          = ',a12,'                  (default)')

140  format ('Boltzmann weighting from min (temp)    = ',f21.8)
145  format ('Boltzmann weighting from min (temp)    = ',f21.8,'         (default)')

120  format ('QNB penalty weight (weight)            = ',f21.8)
125  format ('QNB penalty weight (weight)            = ',f21.8,'         (default)')

end subroutine ffdev_err_nbpnl_ctrl

! ------------------------------------------------------------------------------

end module ffdev_err_nbpnl_control
