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

module ffdev_err_dihedrals_control

use ffdev_constants

contains

! ==============================================================================
! subroutine ffdev_err_dihedrals_ctrl
! ==============================================================================

subroutine ffdev_err_dihedrals_ctrl(fin)

    use ffdev_err_dihedrals_dat
    use ffdev_utils
    use prmfile

    implicit none
    type(PRMFILE_TYPE)  :: fin
    ! --------------------------------------------------------------------------

    write(DEV_OUT,*)
    write(DEV_OUT,10)

    if( .not. prmfile_open_section(fin,'dihedrals') ) then
        write(DEV_OUT,115) prmfile_onoff(EnableDihedralsError)
        write(DEV_OUT,135) prmfile_onoff(PrintDihedralsErrorSummary)
        write(DEV_OUT,125) DihedralsErrorWeight
        write(DEV_OUT,145) prmfile_onoff(OnlyFFOptDihedrals)
        return
    end if

    if( prmfile_get_logical_by_key(fin,'enabled', EnableDihedralsError)) then
        write(DEV_OUT,110) prmfile_onoff(EnableDihedralsError)
    else
        write(DEV_OUT,115) prmfile_onoff(EnableDihedralsError)
    end if
    if( prmfile_get_logical_by_key(fin,'summary', PrintDihedralsErrorSummary)) then
        write(DEV_OUT,130) prmfile_onoff(PrintDihedralsErrorSummary)
    else
        write(DEV_OUT,135) prmfile_onoff(PrintDihedralsErrorSummary)
    end if
    if( prmfile_get_real8_by_key(fin,'weight', DihedralsErrorWeight)) then
        write(DEV_OUT,120) DihedralsErrorWeight
    else
        write(DEV_OUT,125) DihedralsErrorWeight
    end if
    if( prmfile_get_logical_by_key(fin,'onlyffopt', OnlyFFOptDihedrals)) then
        write(DEV_OUT,140) prmfile_onoff(OnlyFFOptDihedrals)
    else
        write(DEV_OUT,145) prmfile_onoff(OnlyFFOptDihedrals)
    end if

 10 format('=== [dihedrals] ================================================================')

110  format ('Dihedrals error (enabled)              = ',a12)
115  format ('Dihedrals error (enabled)              = ',a12,'                  (default)')
130  format ('Print dihedrals error summary (summary)= ',a12)
135  format ('Print dihedrals error summary (summary)= ',a12,'                  (default)')
120  format ('Dihedrals error weight (weight)        = ',f21.8)
125  format ('Dihedrals error weight (weight)        = ',f21.8,'         (default)')
140  format ('Only FFopt dihedrals (onlyffopt)       = ',a12)
145  format ('Only FFopt dihedrals (onlyffopt)       = ',a12,'                  (default)')

end subroutine ffdev_err_dihedrals_ctrl

! ------------------------------------------------------------------------------

end module ffdev_err_dihedrals_control
