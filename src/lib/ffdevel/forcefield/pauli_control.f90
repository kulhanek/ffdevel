! ==============================================================================
! This file is part of FFDevel.
!    Copyright (C) 2019 Petr Kulhanek, kulhanek@chemi.muni.cz
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

module ffdev_pauli_control

use ffdev_sizes
use ffdev_constants

implicit none
contains

!===============================================================================
!-------------------------------------------------------------------------------
!===============================================================================

subroutine ffdev_pauli_ctrl(fin)

    use ffdev_pauli_dat
    use ffdev_pauli
    use prmfile
    use ffdev_utils

    implicit none
    type(PRMFILE_TYPE)  :: fin
    ! --------------------------------------------
    character(80)       :: string
    ! --------------------------------------------------------------------------

    call ffdev_pauli_set_default

    write(DEV_OUT,'(/,a)') '=== [pauli] ===================================================================='

    ! open first section
    if( .not. prmfile_open_section(fin,'pauli') ) then
        write(DEV_OUT,15) prmfile_onoff(pauli_use_numgrid)
        write(DEV_OUT,25) prmfile_onoff(pauli_cache_grid)
        write(DEV_OUT,35) pauli_wf_nsto
        write(DEV_OUT,45) prmfile_onoff(pauli_wf_truncbyn)
        write(DEV_OUT,55) pauli_wf2rho_power
        select case(pauli_wf_form)
            case(PAULI_WF_SV)
                write(DEV_OUT,65) 'SV'
            case(PAULI_WF_STO)
                write(DEV_OUT,65) 'STO'
            case default
                call ffdev_utils_exit(DEV_OUT,1,'Unknown wf_type!')
        end select

        write(DEV_OUT,75) pauli_dens_dpower
        write(DEV_OUT,85) pauli_dens_rpower

        select case(pauli_xfun)
            case(PAULI_XFUN_RHOP)
                write(DEV_OUT,95) 'RHOP'
            case(PAULI_XFUN_KXEG)
                write(DEV_OUT,95) 'KXEG'
            case default
                call ffdev_utils_exit(DEV_OUT,1,'Unknown xfun_type!')
        end select

        write(DEV_OUT,105) pauli_xfun_dpower
        write(DEV_OUT,115) pauli_xfun_kpower
        write(DEV_OUT,125) pauli_xfun_xpower
        write(DEV_OUT,135) pauli_xfun_xfac
        return
    end if

    if( prmfile_get_logical_by_key(fin,'numgrid', pauli_use_numgrid)) then
        write(DEV_OUT,10) prmfile_onoff(pauli_use_numgrid)
    else
        write(DEV_OUT,15) prmfile_onoff(pauli_use_numgrid)
    end if
    if( prmfile_get_logical_by_key(fin,'cache', pauli_cache_grid)) then
        write(DEV_OUT,20) prmfile_onoff(pauli_cache_grid)
    else
        write(DEV_OUT,25) prmfile_onoff(pauli_cache_grid)
    end if

    if( prmfile_get_integer_by_key(fin,'nsto', pauli_wf_nsto)) then
        write(DEV_OUT,30) pauli_wf_nsto
    else
        write(DEV_OUT,35) pauli_wf_nsto
    end if
    if( prmfile_get_logical_by_key(fin,'truncbyn', pauli_wf_truncbyn)) then
        write(DEV_OUT,40) prmfile_onoff(pauli_wf_truncbyn)
    else
        write(DEV_OUT,45) prmfile_onoff(pauli_wf_truncbyn)
    end if
    if( prmfile_get_integer_by_key(fin,'wf2rho_power', pauli_wf2rho_power)) then
        write(DEV_OUT,50) pauli_wf2rho_power
    else
        write(DEV_OUT,55) pauli_wf2rho_power
    end if

    if( prmfile_get_string_by_key(fin,'wf_type', string)) then
        select case(string)
            case('SV')
                pauli_wf_form=PAULI_WF_SV
                write(DEV_OUT,60) trim(string)
            case('STO')
                pauli_wf_form=PAULI_WF_STO
                write(DEV_OUT,60) trim(string)
            case default
                call ffdev_utils_exit(DEV_OUT,1,'Unknown wf_type!')
        end select
    else
        select case(pauli_wf_form)
            case(PAULI_WF_SV)
                write(DEV_OUT,65) 'SV'
            case(PAULI_WF_STO)
                write(DEV_OUT,65) 'STO'
            case default
                call ffdev_utils_exit(DEV_OUT,1,'Unknown wf_type!')
        end select
    end if

    if( prmfile_get_real8_by_key(fin,'pauli_dens_dpower', pauli_dens_dpower)) then
        write(DEV_OUT,70) pauli_dens_dpower
    else
        write(DEV_OUT,75) pauli_dens_dpower
    end if
    if( prmfile_get_real8_by_key(fin,'pauli_dens_rpower', pauli_dens_rpower)) then
        write(DEV_OUT,80) pauli_dens_rpower
    else
        write(DEV_OUT,85) pauli_dens_rpower
    end if

    if( prmfile_get_string_by_key(fin,'xfun_type', string)) then
        select case(string)
            case('SV')
                pauli_xfun=PAULI_XFUN_RHOP
                write(DEV_OUT,90) trim(string)
            case('STO')
                pauli_xfun=PAULI_XFUN_KXEG
                write(DEV_OUT,90) trim(string)
            case default
                call ffdev_utils_exit(DEV_OUT,1,'Unknown xfun_type!')
        end select
    else
        select case(pauli_xfun)
            case(PAULI_XFUN_RHOP)
                write(DEV_OUT,95) 'RHOP'
            case(PAULI_XFUN_KXEG)
                write(DEV_OUT,95) 'KXEG'
            case default
                call ffdev_utils_exit(DEV_OUT,1,'Unknown xfun_type!')
        end select
    end if

    if( prmfile_get_real8_by_key(fin,'pauli_xfun_dpower', pauli_xfun_dpower)) then
        write(DEV_OUT,100) pauli_xfun_dpower
    else
        write(DEV_OUT,105) pauli_xfun_dpower
    end if
    if( prmfile_get_real8_by_key(fin,'pauli_xfun_kpower', pauli_xfun_kpower)) then
        write(DEV_OUT,110) pauli_xfun_kpower
    else
        write(DEV_OUT,115) pauli_xfun_kpower
    end if

    if( prmfile_get_real8_by_key(fin,'pauli_xfun_xpower', pauli_xfun_xpower)) then
        write(DEV_OUT,120) pauli_xfun_xpower
    else
        write(DEV_OUT,125) pauli_xfun_xpower
    end if
    if( prmfile_get_real8_by_key(fin,'pauli_xfun_xfac', pauli_xfun_xfac)) then
        write(DEV_OUT,130) pauli_xfun_xfac
    else
        write(DEV_OUT,135) pauli_xfun_xfac
    end if

    return

 10  format ('Use numgrid for integration (numgrid)  = ',a16)
 15  format ('Use numgrid for integration (numgrid)  = ',a16,'              (default)')
 20  format ('Cache numgrid for integration (cache)  = ',a16)
 25  format ('Cache numgrid for integration (cache)  = ',a16,'              (default)')

 30  format ('Number of STO functions in WF (nsto)   = ',i12)
 35  format ('Number of STO functions in WF (nsto)   = ',i12,'                  (default)')
 40  format ('Truncate STO by n quant num (truncbyn) = ',a16)
 45  format ('Truncate STO by n quant num (truncbyn) = ',a16,'              (default)')
 50  format ('WF to rho transform (wf2rho_power)     = ',i12)
 55  format ('WF to rho transform (wf2rho_power)     = ',i12,'                  (default)')
 60  format ('Type of WF (wf_type)                   = ',a26)
 65  format ('Type of WF (wf_type)                   = ',a26,'    (default)')

 70  format ('pauli_dens_dpower                      = ',f16.8)
 75  format ('pauli_dens_dpower                      = ',f16.8,'              (default)')
 80  format ('pauli_dens_rpower                      = ',f16.8)
 85  format ('pauli_dens_rpower                      = ',f16.8,'              (default)')

 90  format ('Type of XR functional (xfun_type)      = ',a26)
 95  format ('Type of XR functional (xfun_type)      = ',a26,'    (default)')
100  format ('pauli_xfun_dpower                      = ',f16.8)
105  format ('pauli_xfun_dpower                      = ',f16.8,'              (default)')
110  format ('pauli_xfun_kpower                      = ',f16.8)
115  format ('pauli_xfun_kpower                      = ',f16.8,'              (default)')
120  format ('pauli_xfun_xpower                      = ',f16.8)
125  format ('pauli_xfun_xpower                      = ',f16.8,'              (default)')
130  format ('pauli_xfun_xfac                        = ',f16.8)
135  format ('pauli_xfun_xfac                        = ',f16.8,'              (default)')

end subroutine ffdev_pauli_ctrl

!===============================================================================
!-------------------------------------------------------------------------------
!===============================================================================

end module ffdev_pauli_control
