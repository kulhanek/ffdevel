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

module ffdev_err_pbpnl

use ffdev_constants
use ffdev_variables

contains

! ==============================================================================
! subroutine ffdev_err_pbpnl_init
! ==============================================================================

subroutine ffdev_err_pbpnl_init

    use ffdev_err_pbpnl_dat
    use ffdev_errors_dat

    implicit none
    ! --------------------------------------------------------------------------

    EnablePBPnlError        = .false.
    PrintPBPnlErrorSummary  = .false.
    PBPNLMode               = PBPNL_MODE_ALL
    PBPNLSource             = PBPNL_SOURCE_DO
    PBPNLIncludeProbes      = .false.
    PBPnlErrorWeight        = 1.0
    PBPnlErrorWeight1       = 1.0
    PBPnlErrorWeight2       = 0.0
    PBPnlFce                = PBPNL_QUADRATIC

end subroutine ffdev_err_pbpnl_init

! ==============================================================================
! subroutine ffdev_err_pbpnl_error
! ==============================================================================

subroutine ffdev_err_pbpnl_error(error)

    use ffdev_utils
    use ffdev_errors_dat
    use ffdev_err_pbpnl_dat
    use ffdev_parameters_dat
    use ffdev_atomoverlap
    use ffdev_buried_dat
    use ffdev_ip_db

    implicit none
    type(FFERROR_TYPE)  :: error
    ! --------------------------------------------
    integer             :: i
    real(DEVDP)         :: diff,w,pb0
    ! --------------------------------------------------------------------------

    error%pbpnl = 0.0d0

    do i=1,nparams
        if( params(i)%realm .ne. REALM_VDW_PB ) cycle   ! only PB params
        if( params(i)%ti .ne. params(i)%tj ) cycle  ! only like params

        ! weight
        select case(PBPNLMode)
            case(PBPNL_MODE_ALL)
                w = PBPnlErrorWeight1
            case(PBPNL_MODE_BURIED)
                w = PBPnlErrorWeight2 * buried_atoms(params(i)%ti)%weight + &
                    PBPnlErrorWeight1 * (1.0d0 - buried_atoms(params(i)%ti)%weight)
            case(PBPNL_MODE_NOH)
                if( types(params(i)%ti)%z .ne. 1 ) then
                    w = PBPnlErrorWeight1 ! no-H
                else
                    w = PBPnlErrorWeight2 ! H
                end if
            case default
                call ffdev_utils_exit(DEV_ERR,1,'Not implemented mode in ffdev_err_pbpnl_error!')
        end select
        ! probes
        if( PBPNLIncludeProbes .and. types(params(i)%ti)%probe ) then
            w = PBPnlErrorWeight1
        end if

        ! write(*,*) trim(types(params(i)%ti)%name), PBPNLIncludeProbes, types(params(i)%ti)%probe

        select case(PBPNLSource)
            case(PBPNL_SOURCE_DO)
                pb0 = ffdev_atomoverlap_do_bii(params(i)%ti)
            case(PBPNL_SOURCE_IP)
                pb0 = ffdev_bfac_from_ip(params(i)%ti)
            case(PBPNL_SOURCE_IP_XDM)
                pb0 = ffdev_bfac_from_ip_xdm(params(i)%ti)
            case default
                call ffdev_utils_exit(DEV_ERR,1,'Not implemented source in ffdev_err_pbpnl_error!')
        end select

        select case(PBPnlFce)
            case(PBPNL_QUADRATIC)
                diff = params(i)%value - pb0
                error%pbpnl = error%pbpnl + w*diff**2
            case default
                call ffdev_utils_exit(DEV_ERR,1,'Not implemented rst fce in ffdev_err_pbpnl_error!')
        end select
    end do

end subroutine ffdev_err_pbpnl_error

! ==============================================================================
! subroutine ffdev_err_pbpnl_summary
! ==============================================================================

subroutine ffdev_err_pbpnl_summary()

    use ffdev_err_pbpnl_dat
    use ffdev_utils
    use ffdev_parameters_dat
    use ffdev_atomoverlap
    use ffdev_buried_dat
    use ffdev_xdm_dat
    use ffdev_ip_db
    use ffdev_err_pbpnl_control

    implicit none
    integer         :: i
    real(DEVDP)     :: totpnl, pnl, diff, w, pb0
    ! --------------------------------------------------------------------------

    write(DEV_OUT,*)
    write(DEV_OUT,5)
    write(DEV_OUT,7)  trim(ffdev_err_pbpnl_control_mode_to_string(PBPNLMode))
    write(DEV_OUT,8)  trim(ffdev_err_pbpnl_control_source_to_string(PBPNLSource))
    write(DEV_OUT,10)
    write(DEV_OUT,20)

    totpnl = 0.0d0

    do i=1,nparams
        if( params(i)%realm .ne. REALM_VDW_PB ) cycle   ! only PB params
        if( params(i)%ti .ne. params(i)%tj ) cycle  ! only like params

        ! weight
        select case(PBPNLMode)
            case(PBPNL_MODE_ALL)
                w = PBPnlErrorWeight1
            case(PBPNL_MODE_BURIED)
                w = PBPnlErrorWeight2 * buried_atoms(params(i)%ti)%weight + &
                    PBPnlErrorWeight1 * (1.0d0 - buried_atoms(params(i)%ti)%weight)
            case(PBPNL_MODE_NOH)
                if( types(params(i)%ti)%z .eq. 1 ) then
                    w = PBPnlErrorWeight2 ! H
                else
                    w = PBPnlErrorWeight1 ! no-H
                end if
            case default
                call ffdev_utils_exit(DEV_ERR,1,'Not implemented mode in ffdev_err_pbpnl_summary!')
        end select
        ! probes
        if( PBPNLIncludeProbes .and. types(params(i)%ti)%probe ) then
            w = PBPnlErrorWeight1
        end if

        select case(PBPNLSource)
            case(PBPNL_SOURCE_DO)
                pb0 = ffdev_atomoverlap_do_bii(params(i)%ti)
            case(PBPNL_SOURCE_IP)
                pb0 = ffdev_bfac_from_ip(params(i)%ti)
            case(PBPNL_SOURCE_IP_XDM)
                pb0 = ffdev_bfac_from_ip_xdm(params(i)%ti)
            case default
                call ffdev_utils_exit(DEV_ERR,1,'Not implemented source in ffdev_err_pbpnl_error!')
        end select

        select case(PBPnlFce)
            case(PBPNL_QUADRATIC)
                diff = params(i)%value - pb0
                pnl = w*diff**2
            case default
                call ffdev_utils_exit(DEV_ERR,1,'Not implemented in ffdev_err_pbpnl_summary!')
        end select
        totpnl = totpnl + pnl

        write(DEV_OUT,30) i, trim(types(params(i)%ti)%name), pb0, params(i)%value, diff, w, pnl
    end do

    write(DEV_OUT,20)
    write(DEV_OUT,40) totpnl
    write(DEV_OUT,45) totpnl*PBPnlErrorWeight

 5 format('# Pauli Repulsion PB Penalties')
 7 format('# PBPNL Mode: ',A)
 8 format('# PB Source:  ',A)
10 format('# ID Type    PB(Tab)    PB(Opt)       Diff     Weight   Penalty')
20 format('# -- ---- ---------- ---------- ---------- ---------- ----------')
30 format(I4,1X,A4,1X,F10.4,1X,F10.4,1X,F10.5,1X,F10.5,1X,F10.5)
40 format('# Final penalty               =                       ',F10.3)
45 format('# Final penalty (all weights) =                       ',F10.3)

end subroutine ffdev_err_pbpnl_summary

! ------------------------------------------------------------------------------

end module ffdev_err_pbpnl


