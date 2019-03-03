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

module ffdev_err_rmsd_control

use ffdev_constants

contains

! ==============================================================================
! subroutine ffdev_err_rmsd_ctrl
! ==============================================================================

subroutine ffdev_err_rmsd_ctrl(fin)

    use ffdev_err_rmsd_dat
    use ffdev_utils
    use prmfile

    implicit none
    type(PRMFILE_TYPE)  :: fin
    ! --------------------------------------------------------------------------

    write(DEV_OUT,*)
    write(DEV_OUT,10)

    if( .not. prmfile_open_section(fin,'rmsd') ) then
        write(DEV_OUT,115) prmfile_onoff(EnableRMSDError)
        write(DEV_OUT,135) prmfile_onoff(PrintRMSDErrorSummary)
        write(DEV_OUT,125) RMSDErrorWeight
        return
    end if

    if( prmfile_get_logical_by_key(fin,'enabled', EnableRMSDError)) then
        write(DEV_OUT,110) prmfile_onoff(EnableRMSDError)
    else
        write(DEV_OUT,115) prmfile_onoff(EnableRMSDError)
    end if
    if( prmfile_get_logical_by_key(fin,'summary', PrintRMSDErrorSummary)) then
        write(DEV_OUT,130) prmfile_onoff(PrintRMSDErrorSummary)
    else
        write(DEV_OUT,135) prmfile_onoff(PrintRMSDErrorSummary)
    end if
    if( prmfile_get_real8_by_key(fin,'weight', RMSDErrorWeight)) then
        write(DEV_OUT,120) RMSDErrorWeight
    else
        write(DEV_OUT,125) RMSDErrorWeight
    end if

 10 format('=== [rmsd] ==================================================================')

110  format ('RMSD error (enabled)                   = ',a12)
115  format ('RMSD error (enabled)                   = ',a12,'                  (default)')
130  format ('RMSD error summary (summary)           = ',a12)
135  format ('RMSD error summary (summary)           = ',a12,'                  (default)')
120  format ('RMSD error weight (weight)             = ',f21.8)
125  format ('RMSD error weight (weight)             = ',f21.8,'         (default)')

end subroutine ffdev_err_rmsd_ctrl

! ------------------------------------------------------------------------------

end module ffdev_err_rmsd_control
