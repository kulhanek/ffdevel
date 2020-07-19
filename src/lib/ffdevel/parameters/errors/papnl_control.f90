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

module ffdev_err_papnl_control

use ffdev_constants
use ffdev_variables

contains

! ==============================================================================
! subroutine ffdev_err_papnl_ctrl
! ==============================================================================

subroutine ffdev_err_papnl_ctrl(fin)

    use ffdev_err_papnl_dat
    use ffdev_errors_dat
    use ffdev_utils
    use ffdev_errors_utils
    use ffdev_xdm_dat
    use prmfile

    implicit none
    type(PRMFILE_TYPE)          :: fin
    character(80)               :: string
    ! --------------------------------------------------------------------------

    write(DEV_OUT,*)
    write(DEV_OUT,10)

    if( .not. prmfile_open_section(fin,'papnl') ) then
        write(DEV_OUT,115) prmfile_onoff(EnablePAPnlError)
        write(DEV_OUT,135) prmfile_onoff(PrintPAPnlErrorSummary)

        write(DEV_OUT,145) ffdev_err_papnl_control_mode_to_string(PAPNLMode)
        write(DEV_OUT,155) ffdev_err_papnl_control_source_to_string(PAPNLSource)

        write(DEV_OUT,123) PAPnlErrorWeight
        write(DEV_OUT,125) PAPnlErrorWeight1
        write(DEV_OUT,128) PAPnlErrorWeight2
        return
    end if

    if( prmfile_get_logical_by_key(fin,'enabled', EnablePAPnlError)) then
        write(DEV_OUT,110) prmfile_onoff(EnablePAPnlError)
    else
        write(DEV_OUT,115) prmfile_onoff(EnablePAPnlError)
    end if

    if( prmfile_get_logical_by_key(fin,'summary', PrintPAPnlErrorSummary)) then
        write(DEV_OUT,130) prmfile_onoff(PrintPAPnlErrorSummary)
    else
        write(DEV_OUT,135) prmfile_onoff(PrintPAPnlErrorSummary)
    end if

    if( prmfile_get_string_by_key(fin,'mode', string)) then
        PAPNLMode = ffdev_err_papnl_control_mode_from_string(string)
        write(DEV_OUT,140) ffdev_err_papnl_control_mode_to_string(PAPNLMode)
    else
        write(DEV_OUT,145) ffdev_err_papnl_control_mode_to_string(PAPNLMode)
    end if

    if( prmfile_get_string_by_key(fin,'source', string)) then
        PAPNLSource = ffdev_err_papnl_control_source_from_string(string)
        write(DEV_OUT,150) ffdev_err_papnl_control_source_to_string(PAPNLSource)
    else
        write(DEV_OUT,155) ffdev_err_papnl_control_source_to_string(PAPNLSource)
    end if

    if( prmfile_get_real8_by_key(fin,'weight', PAPnlErrorWeight)) then
        write(DEV_OUT,122) PAPnlErrorWeight
    else
        write(DEV_OUT,123) PAPnlErrorWeight
    end if
    if( prmfile_get_real8_by_key(fin,'weight1', PAPnlErrorWeight1)) then
        write(DEV_OUT,120) PAPnlErrorWeight1
    else
        write(DEV_OUT,125) PAPnlErrorWeight1
    end if
    if( prmfile_get_real8_by_key(fin,'weight2', PAPnlErrorWeight2)) then
        write(DEV_OUT,127) PAPnlErrorWeight2
    else
        write(DEV_OUT,128) PAPnlErrorWeight2
    end if

 10 format('=== [papnl] ====================================================================')

110  format ('PA penalty (enabled)                   = ',a12)
115  format ('PA penalty (enabled)                   = ',a12,'                  (default)')
130  format ('PA penalty summary (summary)           = ',a12)
135  format ('PA penalty summary (summary)           = ',a12,'                  (default)')

122  format ('PA penalty weight (weight)             = ',f21.8)
123  format ('PA penalty weight (weight)             = ',f21.8,'         (default)')

120  format ('PA penalty weight#1 (weight1)          = ',f21.8)
125  format ('PA penalty weight#1 (weight1)          = ',f21.8,'         (default)')

127  format ('PA penalty weight#2 (weight2)          = ',f21.8)
128  format ('PA penalty weight#2 (weight2)          = ',f21.8,'         (default)')

140  format ('PA penalty mode (mode)                 = ',a12)
145  format ('PA penalty mode (mode)                 = ',a12,'                  (default)')

150  format ('PA source (source)                     = ',a12)
155  format ('PA source (source)                     = ',a12,'                  (default)')

end subroutine ffdev_err_papnl_ctrl

! ==============================================================================
! subroutine ffdev_err_papnl_control_source_to_string
! ==============================================================================

character(80) function ffdev_err_papnl_control_source_to_string(source)

    use ffdev_utils
    use ffdev_err_papnl_dat

    implicit none
    integer  :: source
    ! --------------------------------------------------------------------------

    select case(source)
        case(PAPNL_SOURCE_DO)
            ffdev_err_papnl_control_source_to_string = 'DO - Density overlaps'
        case default
            call ffdev_utils_exit(DEV_ERR,1,'Not implemented in ffdev_err_papnl_control_source_to_string!')
    end select

end function ffdev_err_papnl_control_source_to_string

! ==============================================================================
! subroutine ffdev_err_papnl_control_source_from_string
! ==============================================================================

integer function ffdev_err_papnl_control_source_from_string(string)

    use ffdev_utils
    use ffdev_err_papnl_dat

    implicit none
    character(*)   :: string
    ! --------------------------------------------------------------------------

    select case(trim(string))
        case('DO')
            ffdev_err_papnl_control_source_from_string = PAPNL_SOURCE_DO
        case default
            call ffdev_utils_exit(DEV_ERR,1,'Not implemented "' // trim(string) &
                                            // '" in ffdev_err_papnl_control_source_from_string!')
    end select

end function ffdev_err_papnl_control_source_from_string

! ==============================================================================
! subroutine ffdev_err_papnl_control_mode_to_string
! ==============================================================================

character(80) function ffdev_err_papnl_control_mode_to_string(mode)

    use ffdev_utils
    use ffdev_err_papnl_dat

    implicit none
    integer  :: mode
    ! --------------------------------------------------------------------------

    select case(mode)
        case(PAPNL_MODE_ALL)
            ffdev_err_papnl_control_mode_to_string = 'ALL'
        case(PAPNL_MODE_BURIED)
            ffdev_err_papnl_control_mode_to_string = 'BURIED'
        case(PAPNL_MODE_NOH)
            ffdev_err_papnl_control_mode_to_string = 'NOH'
        case default
            call ffdev_utils_exit(DEV_ERR,1,'Not implemented in ffdev_err_papnl_control_mode_to_string!')
    end select

end function ffdev_err_papnl_control_mode_to_string

! ==============================================================================
! subroutine ffdev_err_papnl_control_mode_from_string
! ==============================================================================

integer function ffdev_err_papnl_control_mode_from_string(string)

    use ffdev_utils
    use ffdev_err_papnl_dat

    implicit none
    character(*)   :: string
    ! --------------------------------------------------------------------------

    select case(trim(string))
        case('ALL')
            ffdev_err_papnl_control_mode_from_string = PAPNL_MODE_ALL
        case('BURIED')
            ffdev_err_papnl_control_mode_from_string = PAPNL_MODE_BURIED
        case('NOH')
            ffdev_err_papnl_control_mode_from_string = PAPNL_MODE_NOH
        case default
            call ffdev_utils_exit(DEV_ERR,1,'Not implemented "' // trim(string) //'" in ffdev_err_papnl_control_mode_from_string!')
    end select

end function ffdev_err_papnl_control_mode_from_string

! ------------------------------------------------------------------------------

end module ffdev_err_papnl_control
