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

contains

! ==============================================================================
! subroutine ffdev_err_energy_ctrl
! ==============================================================================

subroutine ffdev_err_energy_ctrl(fin)

    use ffdev_err_energy_dat
    use ffdev_errors_dat
    use ffdev_utils
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
        write(DEV_OUT,145) ffdev_err_energy_ctrl_scale_to_string(EnergyErrorMode)
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
        EnergyErrorMode = ffdev_err_energy_ctrl_get_scale_from_string(string)
        write(DEV_OUT,140) ffdev_err_energy_ctrl_scale_to_string(EnergyErrorMode)
    else
        write(DEV_OUT,145) ffdev_err_energy_ctrl_scale_to_string(EnergyErrorMode)
    end if

    errors_calc_ene = EnableEnergyError

 10 format('=== [energy] ===================================================================')

110  format ('Energy error (enabled)                 = ',a12)
115  format ('Energy error (enabled)                 = ',a12,'                  (default)')
130  format ('Print bonds error summary (summary)    = ',a12)
135  format ('Print bonds error summary (summary)    = ',a12,'                  (default)')
120  format ('Energy error weight (weight)           = ',f21.8)
125  format ('Energy error weight (weight)           = ',f21.8,'         (default)')
140  format ('Error scale (scale)                    = ',a12)
145  format ('Error scale (scale)                    = ',a12,'                  (default)')

end subroutine ffdev_err_energy_ctrl

! ==============================================================================
! subroutine ffdev_err_energy_ctrl_scale_to_string
! ==============================================================================

character(80) function ffdev_err_energy_ctrl_scale_to_string(iscale)

    use ffdev_utils
    use ffdev_err_energy_dat

    implicit none
    integer  :: iscale
    ! --------------------------------------------------------------------------

    select case(iscale)
        case(EE_ABS)
            ffdev_err_energy_ctrl_scale_to_string = 'abs (absolute)'
        case(EE_REL)
            ffdev_err_energy_ctrl_scale_to_string = 'rel (relative)'
        case(EE_LOG)
            ffdev_err_energy_ctrl_scale_to_string = 'log (logarithimic)'
        case default
            call ffdev_utils_exit(DEV_OUT,1,'Not implemented in ffdev_err_energy_ctrl_scale_to_string!')
    end select

end function ffdev_err_energy_ctrl_scale_to_string

! ==============================================================================
! function ffdev_err_energy_ctrl_get_scale_from_string
! ==============================================================================

integer function ffdev_err_energy_ctrl_get_scale_from_string(string)

    use ffdev_utils
    use ffdev_err_energy_dat

    implicit none
    character(*)   :: string
    ! --------------------------------------------------------------------------

    select case(trim(string))
        case('abs')
            ffdev_err_energy_ctrl_get_scale_from_string = EE_ABS
        case('rel')
            ffdev_err_energy_ctrl_get_scale_from_string = EE_REL
        case('log')
            ffdev_err_energy_ctrl_get_scale_from_string = EE_LOG
        case default
            call ffdev_utils_exit(DEV_OUT,1,'Not implemented "' // trim(string) // &
                                            '" in ffdev_err_energy_ctrl_get_scale_from_string!')
    end select

end function ffdev_err_energy_ctrl_get_scale_from_string

! ------------------------------------------------------------------------------

end module ffdev_err_energy_control
