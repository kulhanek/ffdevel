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

module ffdev_err_impropers_control

use ffdev_constants
use ffdev_variables

contains

! ==============================================================================
! subroutine ffdev_err_impropers_ctrl
! ==============================================================================

subroutine ffdev_err_impropers_ctrl(fin)

    use ffdev_err_impropers_dat
    use ffdev_utils
    use prmfile

    implicit none
    type(PRMFILE_TYPE)  :: fin
    ! --------------------------------------------------------------------------

    write(DEV_OUT,*)
    write(DEV_OUT,10)

    if( .not. prmfile_open_section(fin,'impropers') ) then
        write(DEV_OUT,115) prmfile_onoff(EnableImpropersError)
        write(DEV_OUT,135) prmfile_onoff(PrintImpropersErrorSummary)
        write(DEV_OUT,125) ImpropersErrorWeight
        write(DEV_OUT,145) prmfile_onoff(ImpropersErrorLockToPhase)
        write(DEV_OUT,155) prmfile_onoff(OnlyFFOptImpropers)
        return
    end if

    if( prmfile_get_logical_by_key(fin,'enabled', EnableImpropersError)) then
        write(DEV_OUT,110) prmfile_onoff(EnableImpropersError)
    else
        write(DEV_OUT,115) prmfile_onoff(EnableImpropersError)
    end if
    if( prmfile_get_logical_by_key(fin,'summary', PrintImpropersErrorSummary)) then
        write(DEV_OUT,130) prmfile_onoff(PrintImpropersErrorSummary)
    else
        write(DEV_OUT,135) prmfile_onoff(PrintImpropersErrorSummary)
    end if
    if( prmfile_get_real8_by_key(fin,'weight', ImpropersErrorWeight)) then
        write(DEV_OUT,120) ImpropersErrorWeight
    else
        write(DEV_OUT,125) ImpropersErrorWeight
    end if
    if( prmfile_get_logical_by_key(fin,'lock2phase', ImpropersErrorLockToPhase)) then
        write(DEV_OUT,140) prmfile_onoff(ImpropersErrorLockToPhase)
    else
        write(DEV_OUT,145) prmfile_onoff(ImpropersErrorLockToPhase)
    end if
    if( prmfile_get_logical_by_key(fin,'onlyffopt', OnlyFFOptImpropers)) then
        write(DEV_OUT,150) prmfile_onoff(OnlyFFOptImpropers)
    else
        write(DEV_OUT,155) prmfile_onoff(OnlyFFOptImpropers)
    end if

 10 format('=== [impropers] ================================================================')

110  format ('Impropers error (enabled)              = ',a12)
115  format ('Impropers error (enabled)              = ',a12,'                  (default)')
130  format ('Print impropers error summary (summary)= ',a12)
135  format ('Print impropers error summary (summary)= ',a12,'                  (default)')
120  format ('Impropers error weight (weight)        = ',f21.8)
125  format ('Impropers error weight (weight)        = ',f21.8,'         (default)')
140  format ('Lock to phase angle (lock2phase)       = ',a12)
145  format ('Lock to phase angle (lock2phase)       = ',a12,'                  (default)')
150  format ('Only FFopt impropers (onlyffopt)       = ',a12)
155  format ('Only FFopt impropers (onlyffopt)       = ',a12,'                  (default)')

end subroutine ffdev_err_impropers_ctrl

! ------------------------------------------------------------------------------

end module ffdev_err_impropers_control
