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

module ffdev_nb2nb_control

use ffdev_constants
use ffdev_variables

contains

! ==============================================================================
! subroutine ffdev_nb2lj_ctrl
! ==============================================================================

subroutine ffdev_nb2lj_ctrl(fin)

    use ffdev_nb2nb_dat
    use ffdev_nb2nb
    use ffdev_utils
    use prmfile

    implicit none
    type(PRMFILE_TYPE)          :: fin
    character(MAX_PATH)         :: buffer
    integer                     :: imode
    ! --------------------------------------------------------------------------

    write(DEV_OUT,*)
    write(DEV_OUT,10)

    if( .not. prmfile_open_section(fin,'nb2lj') ) then
        write(DEV_OUT,115) ffdev_nb2nb_nb2lj_mode_to_string(NB2LJMode)
        write(DEV_OUT,125) NB2LJCutoffR
        write(DEV_OUT,185) NB2LJTemp
        write(DEV_OUT,195) prmfile_onoff(NB2LJWeighted)
        write(DEV_OUT,135) NB2LJIterGS
        write(DEV_OUT,145) NB2LJIterBS
        write(DEV_OUT,155) NB2LJIterOpt
        write(DEV_OUT,165) NB2LJdr
        write(DEV_OUT,205) NB2LJdrPrint
        write(DEV_OUT,215) prmfile_onoff(NB2LJCalcQNBIsoline)
        write(DEV_OUT,225) NB2LJCutoffRQNB
        write(DEV_OUT,175) NBPotPathCore

        call ffdev_nb2nb_initdirs
        return
    end if

    if( prmfile_get_string_by_key(fin,'mode', buffer)) then
        imode = ffdev_nb2nb_nb2lj_mode_from_string(buffer)
        write(DEV_OUT,110) ffdev_nb2nb_nb2lj_mode_to_string(imode)
        NB2LJMode = imode
    else
        write(DEV_OUT,115) ffdev_nb2nb_nb2lj_mode_to_string(NB2LJMode)
    end if

    if( prmfile_get_real8_by_key(fin,'cutoff', NB2LJCutoffR)) then
        write(DEV_OUT,120) NB2LJCutoffR
    else
        write(DEV_OUT,125) NB2LJCutoffR
    end if

    if( prmfile_get_logical_by_key(fin,'weighted', NB2LJWeighted)) then
        write(DEV_OUT,190) prmfile_onoff(NB2LJWeighted)
    else
        write(DEV_OUT,195) prmfile_onoff(NB2LJWeighted)
    end if

    if( prmfile_get_real8_by_key(fin,'temp', NB2LJTemp)) then
        write(DEV_OUT,180) NB2LJTemp
    else
        write(DEV_OUT,185) NB2LJTemp
    end if

    if( prmfile_get_integer_by_key(fin,'itergs', NB2LJIterGS)) then
        write(DEV_OUT,130) NB2LJIterGS
    else
        write(DEV_OUT,135) NB2LJIterGS
    end if

    if( prmfile_get_integer_by_key(fin,'iterbs', NB2LJIterBS)) then
        write(DEV_OUT,140) NB2LJIterBS
    else
        write(DEV_OUT,145) NB2LJIterBS
    end if

    if( prmfile_get_integer_by_key(fin,'iteropt', NB2LJIterOpt)) then
        write(DEV_OUT,150) NB2LJIterOpt
    else
        write(DEV_OUT,155) NB2LJIterOpt
    end if

    if( prmfile_get_real8_by_key(fin,'dr', NB2LJdr)) then
        write(DEV_OUT,160) NB2LJdr
    else
        write(DEV_OUT,165) NB2LJdr
    end if

    if( prmfile_get_real8_by_key(fin,'dr_print', NB2LJdrPrint)) then
        write(DEV_OUT,200) NB2LJdrPrint
    else
        write(DEV_OUT,205) NB2LJdrPrint
    end if

    if( prmfile_get_logical_by_key(fin,'qnb', NB2LJCalcQNBIsoline)) then
        write(DEV_OUT,210) prmfile_onoff(NB2LJCalcQNBIsoline)
    else
        write(DEV_OUT,215) prmfile_onoff(NB2LJCalcQNBIsoline)
    end if

    if( prmfile_get_real8_by_key(fin,'qnbcutoff', NB2LJCutoffRQNB)) then
        write(DEV_OUT,220) NB2LJCutoffRQNB
    else
        write(DEV_OUT,225) NB2LJCutoffRQNB
    end if

    if( prmfile_get_string_by_key(fin,'nbpotpath', NBPotPathCore)) then
        write(DEV_OUT,170) NBPotPathCore
    else
        write(DEV_OUT,175) NBPotPathCore
    end if

    call ffdev_nb2nb_initdirs

 10 format('=== [nb2lj] ====================================================================')

110  format ('NB2LJ mode (mode)                      = ',a12)
115  format ('NB2LJ mode (mode)                      = ',a12,'                  (default)')

120  format ('Cut-off distance (cutoff)              = ',f21.8)
125  format ('Cut-off distance (cutoff)              = ',f21.8,'         (default)')

130  format ('Num of iters in GS (itergs)            = ',i12)
135  format ('Num of iters in GS (itergs)            = ',i12,'                  (default)')

140  format ('Num of iters in BS (iterbs)            = ',i12)
145  format ('Num of iters in BS (iterbs)            = ',i12,'                  (default)')

150  format ('Num of iters in OPT (iteropt)          = ',i12)
155  format ('Num of iters in OPT (iteropt)          = ',i12,'                  (default)')

160  format ('dr for overlay and QNB calc. (dr)      = ',f21.8)
165  format ('dr for overlay and QNB calc. (dr)      = ',f21.8,'         (default)')

200  format ('dr for printing (dr_print)             = ',f21.8)
205  format ('dr for printing (dr_print)             = ',f21.8,'         (default)')

210  format ('Put QNB isoline into NB sum. (qnb)     = ',a12)
215  format ('Put QNB isoline into NB sum. (qnb)     = ',a12,'                  (default)')

220  format ('QNB isoline r-cutoff (qnbcutoff)       = ',f21.8)
225  format ('QNB isoline r-cutoff (qnbcutoff)       = ',f21.8,'         (default)')

170  format ('NB potential summary path (nbpotpath)  = ',a12)
175  format ('NB potential summary path (nbpotpath)  = ',a12,'                  (default)')

180  format ('Temperature factor for weights (temp)  = ',f21.8)
185  format ('Temperature factor for weights (temp)  = ',f21.8,'         (default)')

190  format ('Boltzmann weighting (weighted)         = ',a12)
195  format ('Boltzmann weighting (weighted)         = ',a12,'                  (default)')

end subroutine ffdev_nb2lj_ctrl

! ------------------------------------------------------------------------------

end module ffdev_nb2nb_control
