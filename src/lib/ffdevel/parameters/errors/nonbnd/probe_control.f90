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

module ffdev_err_probe_control

use ffdev_constants
use ffdev_variables

contains

! ==============================================================================
! subroutine ffdev_err_probe_ctrl
! ==============================================================================

subroutine ffdev_err_probe_ctrl(fin)

    use ffdev_err_probe_dat
    use ffdev_errors_dat
    use ffdev_utils
    use ffdev_errors_utils
    use prmfile

    implicit none
    type(PRMFILE_TYPE)          :: fin
    ! --------------------------------------------------------------------------

    write(DEV_OUT,*)
    write(DEV_OUT,10)

    if( .not. prmfile_open_section(fin,'probe') ) then
        write(DEV_OUT,115) prmfile_onoff(EnableProbeError)
        write(DEV_OUT,135) prmfile_onoff(PrintProbeErrorSummary)
        write(DEV_OUT,225) ProbeErrorWeight
        return
    end if

    if( prmfile_get_logical_by_key(fin,'enabled', EnableProbeError)) then
        write(DEV_OUT,110) prmfile_onoff(EnableProbeError)
    else
        write(DEV_OUT,115) prmfile_onoff(EnableProbeError)
    end if
    if( prmfile_get_logical_by_key(fin,'summary', PrintProbeErrorSummary)) then
        write(DEV_OUT,130) prmfile_onoff(PrintProbeErrorSummary)
    else
        write(DEV_OUT,135) prmfile_onoff(PrintProbeErrorSummary)
    end if

    if( prmfile_get_real8_by_key(fin,'weight', ProbeErrorWeight)) then
        write(DEV_OUT,220) ProbeErrorWeight
    else
        write(DEV_OUT,225) ProbeErrorWeight
    end if

 10 format('=== [probe] ====================================================================')

110  format ('Probe error (enabled)                  = ',a12)
115  format ('Probe error (enabled)                  = ',a12,'                  (default)')
130  format ('Print Probe summary (summary)          = ',a12)
135  format ('Print Probe summary (summary)          = ',a12,'                  (default)')

220  format ('Probe error weight (weight)            = ',f21.8)
225  format ('Probe error weight (weight)            = ',f21.8,'         (default)')

end subroutine ffdev_err_probe_ctrl

! ------------------------------------------------------------------------------

end module ffdev_err_probe_control
