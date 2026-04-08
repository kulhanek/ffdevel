! ==============================================================================
! This file is part of FFDevel.
!    Copyright (C) 2026 Petr Kulhanek, kulhanek@chemi.muni.cz
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

module ffdev_err_ihess_control

use ffdev_constants
use ffdev_variables

contains

! ==============================================================================
! subroutine ffdev_err_ihess_ctrl
! ==============================================================================

subroutine ffdev_err_ihess_ctrl(fin)

    use ffdev_err_ihess_dat
    use ffdev_utils
    use prmfile

    implicit none
    type(PRMFILE_TYPE)  :: fin
    ! --------------------------------------------------------------------------

    write(DEV_OUT,*)
    write(DEV_OUT,10)

    if( .not. prmfile_open_section(fin,'ihess') ) then
        write(DEV_OUT,115) prmfile_onoff(EnableIHessError)
        write(DEV_OUT,135) prmfile_onoff(PrintIHessErrorSummary)
        write(DEV_OUT,155) IHessErrorsWeightBonds
        write(DEV_OUT,165) IHessErrorsWeightAngles
        write(DEV_OUT,145) prmfile_onoff(OnlyFFOptIHess)
        write(DEV_OUT,175) prmfile_onoff(IHessErrorsExcludeDihedrals)
        write(DEV_OUT,185) prmfile_onoff(IHessErrorsExcludeImpropers)
        write(DEV_OUT,195) prmfile_onoff(IHessErrorsExcludeNB)
        return
    end if

    if( prmfile_get_logical_by_key(fin,'enabled', EnableIHessError)) then
        write(DEV_OUT,110) prmfile_onoff(EnableIHessError)
    else
        write(DEV_OUT,115) prmfile_onoff(EnableIHessError)
    end if
    if( prmfile_get_logical_by_key(fin,'summary', PrintIHessErrorSummary)) then
        write(DEV_OUT,130) prmfile_onoff(PrintIHessErrorSummary)
    else
        write(DEV_OUT,135) prmfile_onoff(PrintIHessErrorSummary)
    end if
    if( prmfile_get_real8_by_key(fin,'weight_bonds', IHessErrorsWeightBonds)) then
        write(DEV_OUT,150) IHessErrorsWeightBonds
    else
        write(DEV_OUT,155) IHessErrorsWeightBonds
    end if
    if( prmfile_get_real8_by_key(fin,'weight_angles', IHessErrorsWeightAngles)) then
        write(DEV_OUT,160) IHessErrorsWeightAngles
    else
        write(DEV_OUT,165) IHessErrorsWeightAngles
    end if
    if( prmfile_get_logical_by_key(fin,'onlyffopt', OnlyFFOptIHess)) then
        write(DEV_OUT,140) prmfile_onoff(OnlyFFOptIHess)
    else
        write(DEV_OUT,145) prmfile_onoff(OnlyFFOptIHess)
    end if
    if( prmfile_get_logical_by_key(fin,'exclude_dihedrals', IHessErrorsExcludeDihedrals)) then
        write(DEV_OUT,170) prmfile_onoff(IHessErrorsExcludeDihedrals)
    else
        write(DEV_OUT,175) prmfile_onoff(IHessErrorsExcludeDihedrals)
    end if
    if( prmfile_get_logical_by_key(fin,'exclude_impropers', IHessErrorsExcludeImpropers)) then
        write(DEV_OUT,180) prmfile_onoff(IHessErrorsExcludeImpropers)
    else
        write(DEV_OUT,185) prmfile_onoff(IHessErrorsExcludeImpropers)
    end if
    if( prmfile_get_logical_by_key(fin,'exclude_nb', IHessErrorsExcludeNB)) then
        write(DEV_OUT,190) prmfile_onoff(IHessErrorsExcludeNB)
    else
        write(DEV_OUT,195) prmfile_onoff(IHessErrorsExcludeNB)
    end if


 10 format('=== [ihess] ====================================================================')

110  format ('Bonds error (enabled)                  = ',a12)
115  format ('Bonds error (enabled)                  = ',a12,'                  (default)')
130  format ('Print ihess error summary (summary)    = ',a12)
135  format ('Print ihess error summary (summary)    = ',a12,'                  (default)')

150  format ('IHess weight bonds (weight_bonds)      = ',f21.8)
155  format ('IHess weight bonds (weight_bonds)      = ',f21.8,'         (default)')

160  format ('IHess weight angles (weight_angles)    = ',f21.8)
165  format ('IHess weight angles (weight_angles)    = ',f21.8,'         (default)')

140  format ('Only FFopt ihess (onlyffopt)           = ',a12)
145  format ('Only FFopt ihess (onlyffopt)           = ',a12,'                  (default)')

170  format ('Exclude dihedrals (exlude_dihedrals)   = ',a12)
175  format ('Exclude dihedrals (exlude_dihedrals)   = ',a12,'                  (default)')

180  format ('Exclude impropers (exlude_impropers)   = ',a12)
185  format ('Exclude impropers (exlude_impropers)   = ',a12,'                  (default)')

190  format ('Exclude NB (exlude_nb)                 = ',a12)
195  format ('Exclude NB (exlude_nb)                 = ',a12,'                  (default)')

end subroutine ffdev_err_ihess_ctrl

! ------------------------------------------------------------------------------

end module ffdev_err_ihess_control
