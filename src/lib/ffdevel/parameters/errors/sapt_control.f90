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
        write(DEV_OUT,115) prmfile_onoff(EnableSAPTError)
        write(DEV_OUT,135) prmfile_onoff(PrintSAPTErrorSummary)

        write(DEV_OUT,125) SAPTEleErrorWeight
        write(DEV_OUT,225) SAPTRepErrorWeight
        write(DEV_OUT,325) SAPTDispErrorWeight

        write(DEV_OUT,145) prmfile_onoff(SAPTErrorIndToRep)

        errors_calc_sapt = EnableSAPTError
        return
    end if

    if( prmfile_get_logical_by_key(fin,'enabled', EnableSAPTError)) then
        write(DEV_OUT,110) prmfile_onoff(EnableSAPTError)
    else
        write(DEV_OUT,115) prmfile_onoff(EnableSAPTError)
    end if
    if( prmfile_get_logical_by_key(fin,'summary', PrintSAPTErrorSummary)) then
        write(DEV_OUT,130) prmfile_onoff(PrintSAPTErrorSummary)
    else
        write(DEV_OUT,135) prmfile_onoff(PrintSAPTErrorSummary)
    end if

    if( prmfile_get_real8_by_key(fin,'weight_ele', SAPTEleErrorWeight)) then
        write(DEV_OUT,120) SAPTEleErrorWeight
    else
        write(DEV_OUT,125) SAPTEleErrorWeight
    end if
    if( prmfile_get_real8_by_key(fin,'weight_rep', SAPTRepErrorWeight)) then
        write(DEV_OUT,220) SAPTRepErrorWeight
    else
        write(DEV_OUT,225) SAPTRepErrorWeight
    end if
    if( prmfile_get_real8_by_key(fin,'weight_disp', SAPTDispErrorWeight)) then
        write(DEV_OUT,320) SAPTDispErrorWeight
    else
        write(DEV_OUT,325) SAPTDispErrorWeight
    end if
    if( prmfile_get_logical_by_key(fin,'indasrep', SAPTErrorIndToRep)) then
        write(DEV_OUT,140) prmfile_onoff(SAPTErrorIndToRep)
    else
        write(DEV_OUT,145) prmfile_onoff(SAPTErrorIndToRep)
    end if

    errors_calc_sapt = EnableSAPTError

 10 format('=== [sapt] ====================================================================')

110  format ('SAPT error (enabled)                   = ',a12)
115  format ('SAPT error (enabled)                   = ',a12,'                  (default)')
130  format ('Print SAPT summary (summary)           = ',a12)
135  format ('Print SAPT summary (summary)           = ',a12,'                  (default)')

120  format ('SAPT error weight (weight_ele)         = ',f21.8)
125  format ('SAPT error weight (weight_ele)         = ',f21.8,'         (default)')
220  format ('SAPT error weight (weight_rep)         = ',f21.8)
225  format ('SAPT error weight (weight_rep)         = ',f21.8,'         (default)')
320  format ('SAPT error weight (weight_disp)        = ',f21.8)
325  format ('SAPT error weight (weight_disp)        = ',f21.8,'         (default)')

140  format ('Induction as repulsion                 = ',a12)
145  format ('Induction as repulsion                 = ',a12,'                  (default)')

end subroutine ffdev_err_sapt_ctrl

! ------------------------------------------------------------------------------

end module ffdev_err_sapt_control
