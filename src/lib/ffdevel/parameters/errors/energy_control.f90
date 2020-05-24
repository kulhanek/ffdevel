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

module ffdev_err_energy_control

use ffdev_constants
use ffdev_variables

contains

! ==============================================================================
! subroutine ffdev_err_energy_ctrl
! ==============================================================================

subroutine ffdev_err_energy_ctrl(fin)

    use ffdev_err_energy_dat
    use ffdev_errors_dat
    use ffdev_utils
    use ffdev_errors_utils
    use prmfile

    implicit none
    type(PRMFILE_TYPE)          :: fin
    ! --------------------------------------------
    character(PRMFILE_MAX_PATH) :: string
    ! --------------------------------------------------------------------------

    write(DEV_OUT,*)
    write(DEV_OUT,10)

    if( .not. prmfile_open_section(fin,'energy') ) then
        write(DEV_OUT,115) prmfile_onoff(EnableEnergyError)
        write(DEV_OUT,135) prmfile_onoff(PrintEnergyErrorSummary)
        write(DEV_OUT,125) EnergyErrorWeight
        write(DEV_OUT,145) ffdev_errors_utils_scale_to_string(EnergyErrorMode)
        write(DEV_OUT,155) prmfile_onoff(EnableMaxFilter)
        write(DEV_OUT,165) MaxTargetEnergy
        write(DEV_OUT,175) prmfile_onoff(EnableMinFilter)
        write(DEV_OUT,185) MinTargetEnergy
        errors_calc_ene = EnableEnergyError
        return
    end if

    if( prmfile_get_logical_by_key(fin,'enabled', EnableEnergyError)) then
        write(DEV_OUT,110) prmfile_onoff(EnableEnergyError)
    else
        write(DEV_OUT,115) prmfile_onoff(EnableEnergyError)
    end if
    if( prmfile_get_logical_by_key(fin,'summary', PrintEnergyErrorSummary)) then
        write(DEV_OUT,130) prmfile_onoff(PrintEnergyErrorSummary)
    else
        write(DEV_OUT,135) prmfile_onoff(PrintEnergyErrorSummary)
    end if
    if( prmfile_get_real8_by_key(fin,'weight', EnergyErrorWeight)) then
        write(DEV_OUT,120) EnergyErrorWeight
    else
        write(DEV_OUT,125) EnergyErrorWeight
    end if
    if( prmfile_get_string_by_key(fin,'scale', string)) then
        EnergyErrorMode = ffdev_errors_utils_scale_from_string(string)
        write(DEV_OUT,140) ffdev_errors_utils_scale_to_string(EnergyErrorMode)
    else
        write(DEV_OUT,145) ffdev_errors_utils_scale_to_string(EnergyErrorMode)
    end if

    if( prmfile_get_logical_by_key(fin,'maxfilter', EnableMaxFilter)) then
        write(DEV_OUT,150) prmfile_onoff(EnableMaxFilter)
    else
        write(DEV_OUT,155) prmfile_onoff(EnableMaxFilter)
    end if
    if( prmfile_get_real8_by_key(fin,'maxvalue', MaxTargetEnergy)) then
        write(DEV_OUT,160) MaxTargetEnergy
    else
        write(DEV_OUT,165) MaxTargetEnergy
    end if

    if( prmfile_get_logical_by_key(fin,'minfilter', EnableMinFilter)) then
        write(DEV_OUT,170) prmfile_onoff(EnableMinFilter)
    else
        write(DEV_OUT,175) prmfile_onoff(EnableMinFilter)
    end if
    if( prmfile_get_real8_by_key(fin,'minvalue', MinTargetEnergy)) then
        write(DEV_OUT,180) MinTargetEnergy
    else
        write(DEV_OUT,185) MinTargetEnergy
    end if

    errors_calc_ene = EnableEnergyError

 10 format('=== [energy] ===================================================================')

110  format ('Energy error (enabled)                 = ',a12)
115  format ('Energy error (enabled)                 = ',a12,'                  (default)')
130  format ('Print energy error summary (summary)   = ',a12)
135  format ('Print energy error summary (summary)   = ',a12,'                  (default)')
120  format ('Energy error weight (weight)           = ',f21.8)
125  format ('Energy error weight (weight)           = ',f21.8,'         (default)')
140  format ('Error scale (scale)                    = ',a24)
145  format ('Error scale (scale)                    = ',a24,'      (default)')
150  format ('Enable max energy filter (maxfilter)   = ',a12)
155  format ('Enable max energy filter (maxfilter)   = ',a12,'                  (default)')
160  format ('Max target rnergy (maxvalue)           = ',f21.8)
165  format ('Max target rnergy (maxvalue)           = ',f21.8,'         (default)')
170  format ('Enable min energy filter (minfilter)   = ',a12)
175  format ('Enable min energy filter (minfilter)   = ',a12,'                  (default)')
180  format ('Min target rnergy (maxvalue)           = ',f21.8)
185  format ('Min target rnergy (maxvalue)           = ',f21.8,'         (default)')

end subroutine ffdev_err_energy_ctrl

! ------------------------------------------------------------------------------

end module ffdev_err_energy_control
