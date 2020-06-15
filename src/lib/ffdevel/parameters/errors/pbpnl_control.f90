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

module ffdev_err_pbpnl_control

use ffdev_constants
use ffdev_variables

contains

! ==============================================================================
! subroutine ffdev_err_pbpnl_ctrl
! ==============================================================================

subroutine ffdev_err_pbpnl_ctrl(fin)

    use ffdev_err_pbpnl_dat
    use ffdev_errors_dat
    use ffdev_utils
    use ffdev_errors_utils
    use prmfile

    implicit none
    type(PRMFILE_TYPE)          :: fin
    ! --------------------------------------------------------------------------

    write(DEV_OUT,*)
    write(DEV_OUT,10)

    if( .not. prmfile_open_section(fin,'pbpnl') ) then
        write(DEV_OUT,115) prmfile_onoff(EnablePBPnlError)
        write(DEV_OUT,135) prmfile_onoff(PrintPBPnlErrorSummary)

        write(DEV_OUT,125) PBPnlErrorWeight
        write(DEV_OUT,145) prmfile_onoff(PBPNLBuriedAtoms)
        return
    end if

    if( prmfile_get_logical_by_key(fin,'enabled', EnablePBPnlError)) then
        write(DEV_OUT,110) prmfile_onoff(EnablePBPnlError)
    else
        write(DEV_OUT,115) prmfile_onoff(EnablePBPnlError)
    end if

    if( prmfile_get_logical_by_key(fin,'summary', PrintPBPnlErrorSummary)) then
        write(DEV_OUT,130) prmfile_onoff(PrintPBPnlErrorSummary)
    else
        write(DEV_OUT,135) prmfile_onoff(PrintPBPnlErrorSummary)
    end if

    if( prmfile_get_real8_by_key(fin,'weight', PBPnlErrorWeight)) then
        write(DEV_OUT,120) PBPnlErrorWeight
    else
        write(DEV_OUT,125) PBPnlErrorWeight
    end if

    if( prmfile_get_logical_by_key(fin,'buried', PBPNLBuriedAtoms)) then
        write(DEV_OUT,140) prmfile_onoff(PBPNLBuriedAtoms)
    else
        write(DEV_OUT,145) prmfile_onoff(PBPNLBuriedAtoms)
    end if

 10 format('=== [pbpnl] ====================================================================')

110  format ('Charge penalty (enabled)               = ',a12)
115  format ('Charge penalty (enabled)               = ',a12,'                  (default)')
130  format ('Charge penalty summary (summary)       = ',a12)
135  format ('Charge penalty summary (summary)       = ',a12,'                  (default)')

120  format ('Charge penalty weight (weight)         = ',f21.8)
125  format ('Charge penalty weight (weight)         = ',f21.8,'         (default)')

140  format ('Consider buried atoms (buried)         = ',a12)
145  format ('Consider buried atoms (buried)         = ',a12,'                  (default)')

end subroutine ffdev_err_pbpnl_ctrl

! ------------------------------------------------------------------------------

end module ffdev_err_pbpnl_control
