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

module ffdev_err_freqs_control

use ffdev_constants

contains

! ==============================================================================
! subroutine ffdev_err_freqs_ctrl
! ==============================================================================

subroutine ffdev_err_freqs_ctrl(fin)

    use ffdev_err_freqs_dat
    use ffdev_errors_dat
    use ffdev_utils
    use prmfile

    implicit none
    type(PRMFILE_TYPE)  :: fin
    character(80)       :: string
    ! --------------------------------------------------------------------------

    write(DEV_OUT,*)
    write(DEV_OUT,10)

    if( .not. prmfile_open_section(fin,'freqs') ) then
        write(DEV_OUT,115) prmfile_onoff(EnableFreqsError)
        write(DEV_OUT,135) prmfile_onoff(PrintFreqsErrorSummary)
        write(DEV_OUT,145) prmfile_onoff(DebugFreqError)
        write(DEV_OUT,125) FreqsErrorWeight
        write(DEV_OUT,155) FreqsMaxRMSD
        write(DEV_OUT,165) FreqsMaxAngle
        select case(FreqsErrorMode)
            case(FREQS_SUPERIMPOSE_GEO)
                write(DEV_OUT,175) 'bygeo'
            case(FREQS_SUPERIMPOSE_MODE)
                write(DEV_OUT,175) 'bymodes'
        end select
        errors_calc_freq = EnableFreqsError
        return
    end if

    if( prmfile_get_logical_by_key(fin,'enabled', EnableFreqsError)) then
        write(DEV_OUT,110) prmfile_onoff(EnableFreqsError)
    else
        write(DEV_OUT,115) prmfile_onoff(EnableFreqsError)
    end if
    if( prmfile_get_logical_by_key(fin,'summary', PrintFreqsErrorSummary)) then
        write(DEV_OUT,130) prmfile_onoff(PrintFreqsErrorSummary)
    else
        write(DEV_OUT,135) prmfile_onoff(PrintFreqsErrorSummary)
    end if
    if( prmfile_get_logical_by_key(fin,'debug', DebugFreqError)) then
        write(DEV_OUT,140) prmfile_onoff(DebugFreqError)
    else
        write(DEV_OUT,145) prmfile_onoff(DebugFreqError)
    end if
    if( prmfile_get_real8_by_key(fin,'weight', FreqsErrorWeight)) then
        write(DEV_OUT,120) FreqsErrorWeight
    else
        write(DEV_OUT,125) FreqsErrorWeight
    end if
    if( prmfile_get_real8_by_key(fin,'maxrmsd', FreqsMaxRMSD)) then
        write(DEV_OUT,150) FreqsMaxRMSD
    else
        write(DEV_OUT,155) FreqsMaxRMSD
    end if

    if( prmfile_get_real8_by_key(fin,'maxangle', FreqsMaxAngle)) then
        write(DEV_OUT,160) FreqsMaxAngle
    else
        write(DEV_OUT,165) FreqsMaxAngle
    end if

    if( prmfile_get_string_by_key(fin,'mode', string)) then
        select case(string)
            case('bygeo')
                FreqsErrorMode=FREQS_SUPERIMPOSE_GEO
                write(DEV_OUT,170) 'bygeo'
            case('bymodes')
                FreqsErrorMode=FREQS_SUPERIMPOSE_MODE
                write(DEV_OUT,170) 'bymodes'
            case default
                call ffdev_utils_exit(DEV_OUT,1,'Unknown mode specification!')
        end select
    else
        select case(FreqsErrorMode)
            case(FREQS_SUPERIMPOSE_GEO)
                write(DEV_OUT,175) 'bygeo'
            case(FREQS_SUPERIMPOSE_MODE)
                write(DEV_OUT,175) 'bymodes'
        end select
    end if

    errors_calc_freq = EnableFreqsError

 10 format('=== [freqs] ====================================================================')

110  format ('Freq error (enabled)                   = ',a12)
115  format ('Freq error (enabled)                   = ',a12,'                  (default)')
130  format ('Print freq error summary (summary)     = ',a12)
135  format ('Print freq error summary (summary)     = ',a12,'                  (default)')
140  format ('Debug mode (debug)                     = ',a12)
145  format ('Debug mode (debug)                     = ',a12,'                  (default)')
120  format ('Freq error weight (weight)             = ',f21.8)
125  format ('Freq error weight (weight)             = ',f21.8,'         (default)')
150  format ('Max RMSD to trg structure (maxrmsd)    = ',f21.8)
155  format ('Max RMSD to trg structure (maxrmsd)    = ',f21.8,'         (default)')
160  format ('Max angle between nmodes (maxangle)    = ',f21.8)
165  format ('Max angle between nmodes (maxangle)    = ',f21.8,'         (default)')
170  format ('Comparison mode (mode)                 = ',a12)
175  format ('Comparison mode (mode)                 = ',a12,'                  (default)')

end subroutine ffdev_err_freqs_ctrl

! ------------------------------------------------------------------------------

end module ffdev_err_freqs_control
