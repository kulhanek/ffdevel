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

module ffdev_err_pbpnl_control

use ffdev_constants
use ffdev_variables

contains

! ==============================================================================
! subroutine ffdev_err_pbpnl_ctrl
! ==============================================================================

subroutine ffdev_err_pbpnl_ctrl(fin)

    use ffdev_err_pbpnl_dat
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

    if( .not. prmfile_open_section(fin,'pbpnl') ) then
        write(DEV_OUT,115) prmfile_onoff(EnablePBPnlError)
        write(DEV_OUT,135) prmfile_onoff(PrintPBPnlErrorSummary)

        write(DEV_OUT,145) ffdev_err_pbpnl_control_mode_to_string(PBPNLMode)
        write(DEV_OUT,155) ffdev_err_pbpnl_control_source_to_string(PBPNLSource)

        write(DEV_OUT,123) PBPnlErrorWeight
        write(DEV_OUT,125) PBPnlErrorWeight1
        write(DEV_OUT,128) PBPnlErrorWeight2
        return
    end if

    if( prmfile_get_logical_by_key(fin,'enabled', EnablePBPnlError)) then
        write(DEV_OUT,110) prmfile_onoff(EnablePBPnlError)
    else
        write(DEV_OUT,115) prmfile_onoff(EnablePBPnlError)
    end if

    if( prmfile_get_logical_by_key(fin,'summary', PrintPBPnlErrorSummary)) then
        write(DEV_OUT,130) prmfile_onoff(PrintPBPnlErrorSummary)
    else
        write(DEV_OUT,135) prmfile_onoff(PrintPBPnlErrorSummary)
    end if

    if( prmfile_get_string_by_key(fin,'mode', string)) then
        PBPNLMode = ffdev_err_pbpnl_control_mode_from_string(string)
        write(DEV_OUT,140) ffdev_err_pbpnl_control_mode_to_string(PBPNLMode)
    else
        write(DEV_OUT,145) ffdev_err_pbpnl_control_mode_to_string(PBPNLMode)
    end if

    if( prmfile_get_string_by_key(fin,'source', string)) then
        PBPNLSource = ffdev_err_pbpnl_control_source_from_string(string)
        write(DEV_OUT,150) ffdev_err_pbpnl_control_source_to_string(PBPNLSource)
    else
        write(DEV_OUT,155) ffdev_err_pbpnl_control_source_to_string(PBPNLSource)
    end if

    if( prmfile_get_real8_by_key(fin,'weight', PBPnlErrorWeight)) then
        write(DEV_OUT,122) PBPnlErrorWeight
    else
        write(DEV_OUT,123) PBPnlErrorWeight
    end if
    if( prmfile_get_real8_by_key(fin,'weight1', PBPnlErrorWeight1)) then
        write(DEV_OUT,120) PBPnlErrorWeight1
    else
        write(DEV_OUT,125) PBPnlErrorWeight1
    end if
    if( prmfile_get_real8_by_key(fin,'weight2', PBPnlErrorWeight2)) then
        write(DEV_OUT,127) PBPnlErrorWeight2
    else
        write(DEV_OUT,128) PBPnlErrorWeight2
    end if

    if( PBPNLSource .eq. PBPNL_SOURCE_IP_XDM ) then
        if( .not. xdm_data_loaded ) then
            call ffdev_utils_exit(DEV_ERR,1, &
                 'XDM not loaded for PBPNL_SOURCE_IP_XDM in ffdev_err_pbpnl_ctrl!')
        end if
    end if

 10 format('=== [pbpnl] ====================================================================')

110  format ('PB penalty (enabled)                   = ',a12)
115  format ('PB penalty (enabled)                   = ',a12,'                  (default)')
130  format ('PB penalty summary (summary)           = ',a12)
135  format ('PB penalty summary (summary)           = ',a12,'                  (default)')

122  format ('PB penalty weight (weight)             = ',f21.8)
123  format ('PB penalty weight (weight)             = ',f21.8,'         (default)')

120  format ('PB penalty weight#1 (weight1)          = ',f21.8)
125  format ('PB penalty weight#1 (weight1)          = ',f21.8,'         (default)')

127  format ('PB penalty weight#2 (weight2)          = ',f21.8)
128  format ('PB penalty weight#2 (weight2)          = ',f21.8,'         (default)')

140  format ('PB penalty mode (mode)                 = ',a12)
145  format ('PB penalty mode (mode)                 = ',a12,'                  (default)')

150  format ('PB source (source)                     = ',a12)
155  format ('PB source (source)                     = ',a12,'                  (default)')

end subroutine ffdev_err_pbpnl_ctrl

! ==============================================================================
! subroutine ffdev_err_pbpnl_control_source_to_string
! ==============================================================================

character(80) function ffdev_err_pbpnl_control_source_to_string(source)

    use ffdev_utils
    use ffdev_err_pbpnl_dat

    implicit none
    integer  :: source
    ! --------------------------------------------------------------------------

    select case(source)
        case(PBPNL_SOURCE_DO)
            ffdev_err_pbpnl_control_source_to_string = 'DO - Density overlaps'
        case(PBPNL_SOURCE_IP)
            ffdev_err_pbpnl_control_source_to_string = 'IP - Ionization potentials'
        case(PBPNL_SOURCE_IP_XDM)
            ffdev_err_pbpnl_control_source_to_string = 'IP-XDM - Ionization potentials + XDM mods'
        case default
            call ffdev_utils_exit(DEV_ERR,1,'Not implemented in ffdev_err_pbpnl_control_source_to_string!')
    end select

end function ffdev_err_pbpnl_control_source_to_string

! ==============================================================================
! subroutine ffdev_err_pbpnl_control_source_from_string
! ==============================================================================

integer function ffdev_err_pbpnl_control_source_from_string(string)

    use ffdev_utils
    use ffdev_err_pbpnl_dat

    implicit none
    character(*)   :: string
    ! --------------------------------------------------------------------------

    select case(trim(string))
        case('DO')
            ffdev_err_pbpnl_control_source_from_string = PBPNL_SOURCE_DO
        case('IP')
            ffdev_err_pbpnl_control_source_from_string = PBPNL_SOURCE_IP
        case('IP-XDM')
            ffdev_err_pbpnl_control_source_from_string = PBPNL_SOURCE_IP_XDM
        case default
            call ffdev_utils_exit(DEV_ERR,1,'Not implemented "' // trim(string) &
                                            // '" in ffdev_err_pbpnl_control_source_from_string!')
    end select

end function ffdev_err_pbpnl_control_source_from_string

! ==============================================================================
! subroutine ffdev_err_pbpnl_control_mode_to_string
! ==============================================================================

character(80) function ffdev_err_pbpnl_control_mode_to_string(mode)

    use ffdev_utils
    use ffdev_err_pbpnl_dat

    implicit none
    integer  :: mode
    ! --------------------------------------------------------------------------

    select case(mode)
        case(PBPNL_MODE_ALL)
            ffdev_err_pbpnl_control_mode_to_string = 'ALL'
        case(PBPNL_MODE_BURIED)
            ffdev_err_pbpnl_control_mode_to_string = 'BURIED'
        case(PBPNL_MODE_NOH)
            ffdev_err_pbpnl_control_mode_to_string = 'NOH'
        case default
            call ffdev_utils_exit(DEV_ERR,1,'Not implemented in ffdev_err_pbpnl_control_mode_to_string!')
    end select

end function ffdev_err_pbpnl_control_mode_to_string

! ==============================================================================
! subroutine ffdev_err_pbpnl_control_mode_from_string
! ==============================================================================

integer function ffdev_err_pbpnl_control_mode_from_string(string)

    use ffdev_utils
    use ffdev_err_pbpnl_dat

    implicit none
    character(*)   :: string
    ! --------------------------------------------------------------------------

    select case(trim(string))
        case('ALL')
            ffdev_err_pbpnl_control_mode_from_string = PBPNL_MODE_ALL
        case('BURIED')
            ffdev_err_pbpnl_control_mode_from_string = PBPNL_MODE_BURIED
        case('NOH')
            ffdev_err_pbpnl_control_mode_from_string = PBPNL_MODE_NOH
        case default
            call ffdev_utils_exit(DEV_ERR,1,'Not implemented "' // trim(string) //'" in ffdev_err_pbpnl_control_mode_from_string!')
    end select

end function ffdev_err_pbpnl_control_mode_from_string

! ------------------------------------------------------------------------------

end module ffdev_err_pbpnl_control
