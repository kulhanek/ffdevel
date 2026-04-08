! ==============================================================================
! This file is part of FFDevel.
!    Copyright (C) 2024 Petr Kulhanek, kulhanek@chemi.muni.cz
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

module ffdev_err_nbc6_control

use ffdev_constants
use ffdev_variables

contains

! ==============================================================================
! subroutine ffdev_err_nbc6_ctrl
! ==============================================================================

subroutine ffdev_err_nbc6_ctrl(fin)

    use ffdev_err_nbc6_dat
    use ffdev_errors_dat
    use ffdev_utils
    use ffdev_errors_utils
    use prmfile

    implicit none
    type(PRMFILE_TYPE)          :: fin
    ! --------------------------------------------------------------------------

    write(DEV_OUT,*)
    write(DEV_OUT,10)

    if( .not. prmfile_open_section(fin,'nbc6') ) then
        write(DEV_OUT,115) prmfile_onoff(EnableNBC6Error)
        write(DEV_OUT,135) prmfile_onoff(PrintNBC6ErrorSummary)
        write(DEV_OUT,125) NBC6ErrorWeight
        write(DEV_OUT,145) prmfile_onoff(NBC6BurriedOnly)
        return
    end if

    if( prmfile_get_logical_by_key(fin,'enabled', EnableNBC6Error)) then
        write(DEV_OUT,110) prmfile_onoff(EnableNBC6Error)
    else
        write(DEV_OUT,115) prmfile_onoff(EnableNBC6Error)
    end if

    if( prmfile_get_logical_by_key(fin,'summary', PrintNBC6ErrorSummary)) then
        write(DEV_OUT,130) prmfile_onoff(PrintNBC6ErrorSummary)
    else
        write(DEV_OUT,135) prmfile_onoff(PrintNBC6ErrorSummary)
    end if

    if( prmfile_get_real8_by_key(fin,'weight', NBC6ErrorWeight)) then
        write(DEV_OUT,120) NBC6ErrorWeight
    else
        write(DEV_OUT,125) NBC6ErrorWeight
    end if

    if( prmfile_get_logical_by_key(fin,'buried', NBC6BurriedOnly)) then
        write(DEV_OUT,140) prmfile_onoff(NBC6BurriedOnly)
    else
        write(DEV_OUT,145) prmfile_onoff(NBC6BurriedOnly)
    end if

    if( prmfile_get_logical_by_key(fin,'sqrt16', NBC6Sqrt16)) then
        write(DEV_OUT,150) prmfile_onoff(NBC6Sqrt16)
    else
        write(DEV_OUT,155) prmfile_onoff(NBC6Sqrt16)
    end if

    if( prmfile_get_integer_by_key(fin,'C6eff', NBC6Eff)) then
        write(DEV_OUT,160) NBC6Eff
    else
        write(DEV_OUT,165) NBC6Eff
    end if

 10 format('=== [nbc6] =====================================================================')

110  format ('LJ C6 penalty (enabled)                = ',a12)
115  format ('LJ C6 penalty (enabled)                = ',a12,'                  (default)')
130  format ('LJ C6 penalty summary (summary)        = ',a12)
135  format ('LJ C6 penalty summary (summary)        = ',a12,'                  (default)')

120  format ('LJ C6 penalty weight (weight)          = ',f21.8)
125  format ('LJ C6 penalty weight (weight)          = ',f21.8,'         (default)')

140  format ('LJ C6 penalty only buried (buried)     = ',a12)
145  format ('LJ C6 penalty only buried (buried)     = ',a12,'                  (default)')

150  format ('LJ C6 use err**(1/6) as error (sqrt16) = ',a12)
155  format ('LJ C6 use err**(1/6) as error (sqrt16) = ',a12,'                  (default)')

160  format ('Use C6eff calc from C6,C8,C10 ((C6eff) = ',i12)
165  format ('Use C6eff calc from C6,C8,C10 ((C6eff) = ',i12,'                  (default)')

end subroutine ffdev_err_nbc6_ctrl

! ------------------------------------------------------------------------------

end module ffdev_err_nbc6_control
