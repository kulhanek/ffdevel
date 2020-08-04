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

module ffdev_err_papnl

use ffdev_constants
use ffdev_variables

contains

! ==============================================================================
! subroutine ffdev_err_papnl_init
! ==============================================================================

subroutine ffdev_err_papnl_init

    use ffdev_err_papnl_dat
    use ffdev_errors_dat

    implicit none
    ! --------------------------------------------------------------------------

    EnablePAPnlError        = .false.
    PrintPAPnlErrorSummary  = .false.
    PAPNLMode               = PAPNL_MODE_ALL
    PAPNLSource             = PAPNL_SOURCE_DO
    PAPnlErrorWeight        = 1.0
    PAPnlErrorWeight1       = 1.0
    PAPnlErrorWeight2       = 0.0
    PAPnlFce                = PAPNL_QUADRATIC

end subroutine ffdev_err_papnl_init

! ==============================================================================
! subroutine ffdev_err_papnl_error
! ==============================================================================

subroutine ffdev_err_papnl_error(error)

    use ffdev_utils
    use ffdev_errors_dat
    use ffdev_err_papnl_dat
    use ffdev_parameters_dat
    use ffdev_atomicdata
    use ffdev_buried_dat
    use ffdev_ip_db

    implicit none
    type(FFERROR_TYPE)  :: error
    ! --------------------------------------------
    integer             :: i
    real(DEVDP)         :: diff,w,pa0
    ! --------------------------------------------------------------------------

    error%papnl = 0.0d0

    do i=1,nparams
        if( params(i)%realm .ne. REALM_VDW_PA ) cycle   ! only PB params
        if( params(i)%ti .ne. params(i)%tj ) cycle  ! only like params

        ! weight
        select case(PAPNLMode)
            case(PAPNL_MODE_ALL)
                w = PAPnlErrorWeight1
            case(PAPNL_MODE_BURIED)
                w = PAPnlErrorWeight2 * buried_atoms(params(i)%ti)%weight + &
                    PAPnlErrorWeight1 * (1.0d0 - buried_atoms(params(i)%ti)%weight)
            case(PAPNL_MODE_NOH)
                if( types(params(i)%ti)%z .ne. 1 ) then
                    w = PAPnlErrorWeight1 ! no-H
                else
                    w = PAPnlErrorWeight2 ! H
                end if
            case default
                call ffdev_utils_exit(DEV_ERR,1,'Not implemented mode in ffdev_err_papnl_error!')
        end select

        select case(PAPNLSource)
            case(PAPNL_SOURCE_DO)
                pa0 = ffdev_atomicdata_do_aii(params(i)%ti) + pauli_k
            case default
                call ffdev_utils_exit(DEV_ERR,1,'Not implemented source in ffdev_err_papnl_error!')
        end select

        select case(PAPnlFce)
            case(PAPNL_QUADRATIC)
                diff = params(i)%value - pa0
                error%papnl = error%papnl + w*diff**2
            case default
                call ffdev_utils_exit(DEV_ERR,1,'Not implemented rst fce in ffdev_err_papnl_error!')
        end select
    end do

end subroutine ffdev_err_papnl_error

! ==============================================================================
! subroutine ffdev_err_papnl_summary
! ==============================================================================

subroutine ffdev_err_papnl_summary()

    use ffdev_err_papnl_dat
    use ffdev_utils
    use ffdev_parameters_dat
    use ffdev_atomicdata
    use ffdev_buried_dat
    use ffdev_xdm_dat
    use ffdev_ip_db
    use ffdev_err_papnl_control

    implicit none
    integer         :: i
    real(DEVDP)     :: totpnl, pnl, diff, w, pa0
    ! --------------------------------------------------------------------------

    write(DEV_OUT,*)
    write(DEV_OUT,5)
    write(DEV_OUT,7)  trim(ffdev_err_papnl_control_mode_to_string(PAPNLMode))
    write(DEV_OUT,8)  trim(ffdev_err_papnl_control_source_to_string(PAPNLSource))
    write(DEV_OUT,10)
    write(DEV_OUT,20)

    totpnl = 0.0d0

    do i=1,nparams
        if( params(i)%realm .ne. REALM_VDW_PA ) cycle   ! only PB params
        if( params(i)%ti .ne. params(i)%tj ) cycle  ! only like params

        ! weight
        select case(PAPNLMode)
            case(PAPNL_MODE_ALL)
                w = PAPnlErrorWeight1
            case(PAPNL_MODE_BURIED)
                w = PAPnlErrorWeight2 * buried_atoms(params(i)%ti)%weight + &
                    PAPnlErrorWeight1 * (1.0d0 - buried_atoms(params(i)%ti)%weight)
            case(PAPNL_MODE_NOH)
                if( types(params(i)%ti)%z .eq. 1 ) then
                    w = PAPnlErrorWeight2 ! H
                else
                    w = PAPnlErrorWeight1 ! no-H
                end if
            case default
                call ffdev_utils_exit(DEV_ERR,1,'Not implemented mode in ffdev_err_papnl_summary!')
        end select

        select case(PAPNLSource)
            case(PAPNL_SOURCE_DO)
                pa0 = ffdev_atomicdata_do_aii(params(i)%ti) + pauli_k
            case default
                call ffdev_utils_exit(DEV_ERR,1,'Not implemented source in ffdev_err_papnl_error!')
        end select

        select case(PAPnlFce)
            case(PAPNL_QUADRATIC)
                diff = params(i)%value - pa0
                pnl = w*diff**2
            case default
                call ffdev_utils_exit(DEV_ERR,1,'Not implemented in ffdev_err_papnl_summary!')
        end select
        totpnl = totpnl + pnl

        write(DEV_OUT,30) i, trim(types(params(i)%ti)%name), pa0, params(i)%value, diff, w, pnl
    end do

    write(DEV_OUT,20)
    write(DEV_OUT,40) totpnl
    write(DEV_OUT,45) totpnl*PAPnlErrorWeight

 5 format('# Pauli Repulsion PA Penalties')
 7 format('# PAPNL Mode: ',A)
 8 format('# PA Source:  ',A)
10 format('# ID Type    PA(Tab)    PA(Opt)       Diff     Weight   Penalty')
20 format('# -- ---- ---------- ---------- ---------- ---------- ----------')
30 format(I4,1X,A4,1X,F10.4,1X,F10.4,1X,F10.5,1X,F10.5,1X,F10.5,1X,F10.5)
40 format('# Final penalty               =                       ',F10.3)
45 format('# Final penalty (all weights) =                       ',F10.3)

end subroutine ffdev_err_papnl_summary

! ------------------------------------------------------------------------------

end module ffdev_err_papnl


