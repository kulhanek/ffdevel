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
! subroutine ffdev_nb2nb_ctrl
! ==============================================================================

subroutine ffdev_nb2nb_ctrl(fin)

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

    if( .not. prmfile_open_section(fin,'nb2nb') ) then
        write(DEV_OUT,115) ffdev_nb2nb_nb2nb_mode_to_string(NB2NBMode)
        write(DEV_OUT,125) NB2NBCutoffR
        write(DEV_OUT,185) NB2NBTemp
        write(DEV_OUT,195) prmfile_onoff(NB2NBWeighted)
        write(DEV_OUT,135) NB2NBIterGS
        write(DEV_OUT,145) NB2NBIterBS
        write(DEV_OUT,155) NB2NBIterOpt
        write(DEV_OUT,165) NB2NBdr
        write(DEV_OUT,205) NB2NBdrPrint
        write(DEV_OUT,215) prmfile_onoff(NB2NBCalcQNBIsoline)
        write(DEV_OUT,175) NBPotPathCore

        call ffdev_nb2nb_initdirs
        return
    end if

    if( prmfile_get_string_by_key(fin,'mode', buffer)) then
        imode = ffdev_nb2nb_nb2nb_mode_from_string(buffer)
        write(DEV_OUT,110) ffdev_nb2nb_nb2nb_mode_to_string(imode)
        NB2NBMode = imode
    else
        write(DEV_OUT,115) ffdev_nb2nb_nb2nb_mode_to_string(NB2NBMode)
    end if

    if( prmfile_get_real8_by_key(fin,'cutoff', NB2NBCutoffR)) then
        write(DEV_OUT,120) NB2NBCutoffR
    else
        write(DEV_OUT,125) NB2NBCutoffR
    end if

    if( prmfile_get_logical_by_key(fin,'weighted', NB2NBWeighted)) then
        write(DEV_OUT,190) prmfile_onoff(NB2NBWeighted)
    else
        write(DEV_OUT,195) prmfile_onoff(NB2NBWeighted)
    end if

    if( prmfile_get_real8_by_key(fin,'temp', NB2NBTemp)) then
        write(DEV_OUT,180) NB2NBTemp
    else
        write(DEV_OUT,185) NB2NBTemp
    end if

    if( prmfile_get_integer_by_key(fin,'itergs', NB2NBIterGS)) then
        write(DEV_OUT,130) NB2NBIterGS
    else
        write(DEV_OUT,135) NB2NBIterGS
    end if

    if( prmfile_get_integer_by_key(fin,'iterbs', NB2NBIterBS)) then
        write(DEV_OUT,140) NB2NBIterBS
    else
        write(DEV_OUT,145) NB2NBIterBS
    end if

    if( prmfile_get_integer_by_key(fin,'iteropt', NB2NBIterOpt)) then
        write(DEV_OUT,150) NB2NBIterOpt
    else
        write(DEV_OUT,155) NB2NBIterOpt
    end if

    if( prmfile_get_real8_by_key(fin,'dr', NB2NBdr)) then
        write(DEV_OUT,160) NB2NBdr
    else
        write(DEV_OUT,165) NB2NBdr
    end if

    if( prmfile_get_real8_by_key(fin,'dr_print', NB2NBdrPrint)) then
        write(DEV_OUT,200) NB2NBdrPrint
    else
        write(DEV_OUT,205) NB2NBdrPrint
    end if

    if( prmfile_get_logical_by_key(fin,'qnb', NB2NBCalcQNBIsoline)) then
        write(DEV_OUT,210) prmfile_onoff(NB2NBCalcQNBIsoline)
    else
        write(DEV_OUT,215) prmfile_onoff(NB2NBCalcQNBIsoline)
    end if

    if( prmfile_get_string_by_key(fin,'nbpotpath', NBPotPathCore)) then
        write(DEV_OUT,170) NBPotPathCore
    else
        write(DEV_OUT,175) NBPotPathCore
    end if

    call ffdev_nb2nb_initdirs
    call ffdev_nb2nb_setup_gaussian_quadrature

 10 format('=== [nb2nb] ====================================================================')

110  format ('NB2NB mode (mode)                      = ',a12)
115  format ('NB2NB mode (mode)                      = ',a12,'                  (default)')

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

170  format ('NB potential summary path (nbpotpath)  = ',a12)
175  format ('NB potential summary path (nbpotpath)  = ',a12,'                  (default)')

180  format ('Temperature factor for weights (temp)  = ',f21.8)
185  format ('Temperature factor for weights (temp)  = ',f21.8,'         (default)')

190  format ('Boltzmann weighting (weighted)         = ',a12)
195  format ('Boltzmann weighting (weighted)         = ',a12,'                  (default)')

end subroutine ffdev_nb2nb_ctrl

! ------------------------------------------------------------------------------

end module ffdev_nb2nb_control
