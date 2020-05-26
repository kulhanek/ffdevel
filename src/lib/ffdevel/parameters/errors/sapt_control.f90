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

module ffdev_err_sapt_control

use ffdev_constants
use ffdev_variables

contains

! ==============================================================================
! subroutine ffdev_err_sapt_ctrl
! ==============================================================================

subroutine ffdev_err_sapt_ctrl(fin)

    use ffdev_err_sapt_dat
    use ffdev_errors_dat
    use ffdev_utils
    use ffdev_errors_utils
    use prmfile

    implicit none
    type(PRMFILE_TYPE)          :: fin
    ! --------------------------------------------------------------------------

    write(DEV_OUT,*)
    write(DEV_OUT,10)

    if( .not. prmfile_open_section(fin,'sapt') ) then
        write(DEV_OUT,115) prmfile_onoff(EnableSAPT0Error)
        write(DEV_OUT,135) prmfile_onoff(PrintSAPT0ErrorSummary)

        write(DEV_OUT,125) SAPT0EleErrorWeight
        write(DEV_OUT,225) SAPT0RepErrorWeight
        write(DEV_OUT,325) SAPT0DispErrorWeight

        write(DEV_OUT,145) prmfile_onoff(SAPT0ErrorIndToRep)

        errors_calc_sapt = EnableSAPT0Error
        return
    end if

    if( prmfile_get_logical_by_key(fin,'enabled', EnableSAPT0Error)) then
        write(DEV_OUT,110) prmfile_onoff(EnableSAPT0Error)
    else
        write(DEV_OUT,115) prmfile_onoff(EnableSAPT0Error)
    end if
    if( prmfile_get_logical_by_key(fin,'summary', PrintSAPT0ErrorSummary)) then
        write(DEV_OUT,130) prmfile_onoff(PrintSAPT0ErrorSummary)
    else
        write(DEV_OUT,135) prmfile_onoff(PrintSAPT0ErrorSummary)
    end if

    if( prmfile_get_real8_by_key(fin,'weight_ele', SAPT0EleErrorWeight)) then
        write(DEV_OUT,120) SAPT0EleErrorWeight
    else
        write(DEV_OUT,125) SAPT0EleErrorWeight
    end if
    if( prmfile_get_real8_by_key(fin,'weight_rep', SAPT0RepErrorWeight)) then
        write(DEV_OUT,220) SAPT0RepErrorWeight
    else
        write(DEV_OUT,225) SAPT0RepErrorWeight
    end if
    if( prmfile_get_real8_by_key(fin,'weight_disp', SAPT0DispErrorWeight)) then
        write(DEV_OUT,320) SAPT0DispErrorWeight
    else
        write(DEV_OUT,325) SAPT0DispErrorWeight
    end if
    if( prmfile_get_logical_by_key(fin,'indasrep', SAPT0ErrorIndToRep)) then
        write(DEV_OUT,140) prmfile_onoff(SAPT0ErrorIndToRep)
    else
        write(DEV_OUT,145) prmfile_onoff(SAPT0ErrorIndToRep)
    end if

    errors_calc_sapt = EnableSAPT0Error

 10 format('=== [sapt] ===================================================================')

110  format ('SAPT0 error (enabled)                  = ',a12)
115  format ('SAPT0 error (enabled)                  = ',a12,'                  (default)')
130  format ('Print SAPT0 summary (summary)          = ',a12)
135  format ('Print SAPT0 summary (summary)          = ',a12,'                  (default)')

120  format ('SAPT0 error weight (weight_ele)        = ',f21.8)
125  format ('SAPT0 error weight (weight_ele)        = ',f21.8,'         (default)')
220  format ('SAPT0 error weight (weight_rep)        = ',f21.8)
225  format ('SAPT0 error weight (weight_rep)        = ',f21.8,'         (default)')
320  format ('SAPT0 error weight (weight_disp)       = ',f21.8)
325  format ('SAPT0 error weight (weight_disp)       = ',f21.8,'         (default)')

140  format ('Induction as repulsion                 = ',a12)
145  format ('Induction as repulsion                 = ',a12,'                  (default)')

end subroutine ffdev_err_sapt_ctrl

! ------------------------------------------------------------------------------

end module ffdev_err_sapt_control
