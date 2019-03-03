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

module ffdev_err_nbdists_control

use ffdev_constants

contains

! ==============================================================================
! subroutine ffdev_err_nbdists_ctrl
! ==============================================================================

subroutine ffdev_err_nbdists_ctrl(fin)

    use ffdev_err_nbdists_dat
    use ffdev_utils
    use prmfile

    implicit none
    type(PRMFILE_TYPE)  :: fin
    ! --------------------------------------------------------------------------

    write(DEV_OUT,*)
    write(DEV_OUT,10)

    if( .not. prmfile_open_section(fin,'nbdists') ) then
        write(DEV_OUT,115) prmfile_onoff(EnableNBDistsError)
        write(DEV_OUT,135) prmfile_onoff(PrintNBDistsErrorSummary)
        write(DEV_OUT,125) NBDistsErrorWeight
        write(DEV_OUT,65) NBDistanceSWPosition
        write(DEV_OUT,75) NBDistanceSWAlpha
        return
    end if

    if( prmfile_get_logical_by_key(fin,'enabled', EnableNBDistsError)) then
        write(DEV_OUT,110) prmfile_onoff(EnableNBDistsError)
    else
        write(DEV_OUT,115) prmfile_onoff(EnableNBDistsError)
    end if
    if( prmfile_get_logical_by_key(fin,'summary', PrintNBDistsErrorSummary)) then
        write(DEV_OUT,130) prmfile_onoff(PrintNBDistsErrorSummary)
    else
        write(DEV_OUT,135) prmfile_onoff(PrintNBDistsErrorSummary)
    end if
    if( prmfile_get_real8_by_key(fin,'weight', NBDistsErrorWeight)) then
        write(DEV_OUT,120) NBDistsErrorWeight
    else
        write(DEV_OUT,125) NBDistsErrorWeight
    end if

    if( prmfile_get_real8_by_key(fin,'swr0', NBDistanceSWPosition)) then
        write(DEV_OUT,60) NBDistanceSWPosition
    else
        write(DEV_OUT,65) NBDistanceSWPosition
    end if

    if( prmfile_get_real8_by_key(fin,'swa', NBDistanceSWAlpha)) then
        write(DEV_OUT,60) NBDistanceSWAlpha
    else
        write(DEV_OUT,65) NBDistanceSWAlpha
    end if

 10 format('=== [nbdists] ==================================================================')

110  format ('NB distance error (enabled)            = ',a12)
115  format ('NB distance error (enabled)            = ',a12,'                  (default)')
130  format ('Print NB dist error summary (summary)  = ',a12)
135  format ('Print NB dist error summary (summary)  = ',a12,'                  (default)')

120  format ('NB distance error weight (weight)      = ',f21.8)
125  format ('NB distance error weight (weight)      = ',f21.8,'         (default)')

 60  format ('NB distance switch r0 (swr0)           = ',f12.6)
 65  format ('NB distance switch r0 (swr0)           = ',f12.6,'                  (default)')

 70  format ('NB distance switch alpha (swa)         = ',f12.6)
 75  format ('NB distance switch alpha (swa)         = ',f12.6,'                  (default)')

end subroutine ffdev_err_nbdists_ctrl

! ------------------------------------------------------------------------------

end module ffdev_err_nbdists_control
