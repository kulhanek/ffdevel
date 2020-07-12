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
    PBPnlErrorWeight        = 1.0
    PBPnlFce                = PB_PNL_QUADRATIC
    PBPNLBuriedAtoms        = .false.

end subroutine ffdev_err_pbpnl_init

! ==============================================================================
! subroutine ffdev_err_pbpnl_error
! ==============================================================================

subroutine ffdev_err_pbpnl_error(error)

    use ffdev_utils
    use ffdev_errors_dat
    use ffdev_err_pbpnl_dat
    use ffdev_parameters_dat
    use ffdev_densoverlap
    use ffdev_buried_dat

    implicit none
    type(FFERROR_TYPE)  :: error
    ! --------------------------------------------
    integer             :: i
    real(DEVDP)         :: diff,w
    ! --------------------------------------------------------------------------

    error%pbpnl = 0.0d0

    do i=1,nparams
        if( params(i)%realm .ne. REALM_VDW_PB ) cycle   ! only PB params
        if( params(i)%ti .ne. params(i)%tj ) cycle  ! only like params
        w = 1.0d0
        if( PBPNLBuriedAtoms ) then
            w = 1.0d0 - buried_atoms(params(i)%ti)%weight
        end if
        select case(PBPnlFce)
            case(PB_PNL_QUADRATIC)
                diff = params(i)%value - ffdev_densoverlap_b(params(i)%ti)
                error%pbpnl = error%pbpnl + w*diff**2
            case default
                call ffdev_utils_exit(DEV_ERR,1,'Not implemented in ffdev_err_pbpnl_error!')
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
    use ffdev_densoverlap
    use ffdev_buried_dat

    implicit none
    logical         :: printsum
    integer         :: i
    real(DEVDP)     :: totpnl, pnl, diff, w
    ! --------------------------------------------------------------------------

    printsum = .false.
    do i=1,nparams
        if( params(i)%realm .ne. REALM_VDW_PB ) cycle   ! only PB params
        if( params(i)%ti .ne. params(i)%tj ) cycle  ! only like params
        printsum = .true.
    end do
    if( .not. printsum ) return

    write(DEV_OUT,*)
    write(DEV_OUT,5)
    write(DEV_OUT,10)
    write(DEV_OUT,20)

    totpnl = 0.0d0

    do i=1,nparams
        if( params(i)%realm .ne. REALM_VDW_PB ) cycle   ! only PB params
        if( params(i)%ti .ne. params(i)%tj ) cycle  ! only like params

        w = 1.0d0
        if( PBPNLBuriedAtoms ) then
            w = 1.0d0 - buried_atoms(params(i)%ti)%weight
        end if

        select case(PBPnlFce)
            case(PB_PNL_QUADRATIC)
                diff = params(i)%value - ffdev_densoverlap_b(params(i)%ti)
                pnl = w*diff**2
            case default
                call ffdev_utils_exit(DEV_ERR,1,'Not implemented in ffdev_err_pbpnl_summary!')
        end select
        totpnl = totpnl + pnl

        write(DEV_OUT,30) i, trim(types(params(i)%ti)%name), ffdev_densoverlap_b(params(i)%ti), params(i)%value, diff, pnl
    end do

    write(DEV_OUT,20)
    write(DEV_OUT,40) totpnl
    write(DEV_OUT,45) totpnl*PBPnlErrorWeight

 5 format('# Pauli Repulsion PB Penalties')
10 format('# ID Type    PB(Tab)    PB(Opt)       Diff    Penalty')
20 format('# -- ---- ---------- ---------- ---------- ----------')
30 format(I4,1X,A4,1X,F10.4,1X,F10.4,1X,F10.4,1X,F10.4)
40 format('# Final penalty               =            ',F10.3)
45 format('# Final penalty (all weights) =            ',F10.3)

end subroutine ffdev_err_pbpnl_summary

! ------------------------------------------------------------------------------

end module ffdev_err_pbpnl


