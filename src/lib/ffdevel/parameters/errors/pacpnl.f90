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

module ffdev_err_pacpnl

use ffdev_constants
use ffdev_variables

contains

! ==============================================================================
! subroutine ffdev_err_pacpnl_init
! ==============================================================================

subroutine ffdev_err_pacpnl_init

    use ffdev_err_pacpnl_dat
    use ffdev_errors_dat

    implicit none
    ! --------------------------------------------------------------------------

    EnablePACPnlError        = .false.
    PrintPACPnlErrorSummary  = .false.
    PACPnlErrorWeight        = 1.0
    PACPnlFce                = CHRG_PNL_QUADRATIC

end subroutine ffdev_err_pacpnl_init

! ==============================================================================
! subroutine ffdev_err_pacpnl_error
! ==============================================================================

subroutine ffdev_err_pacpnl_error(error)

    use ffdev_targetset
    use ffdev_targetset_dat
    use ffdev_utils
    use ffdev_geometry
    use ffdev_errors_dat
    use ffdev_err_pacpnl_dat

    implicit none
    type(FFERROR_TYPE)  :: error
    ! --------------------------------------------
    integer             :: i,j,k
    real(DEVDP)         :: diff
    ! --------------------------------------------------------------------------

    error%pacpnl = 0.0d0

    do i=1,nsets
        do j=1,sets(i)%ngeos
            if( .not. sets(i)%geo(j)%sup_chrg_loaded ) cycle
            do k=1,sets(i)%top%natoms
                select case(PACPnlFce)
                    case(CHRG_PNL_QUADRATIC)
                        diff = sets(i)%geo(j)%sup_chrg(k) - sets(i)%top%atoms(k)%charge
                        error%pacpnl = error%pacpnl + diff**2
                    case default
                        call ffdev_utils_exit(DEV_ERR,1,'Not implemented in ffdev_err_pacpnl_error!')
                end select
            end do
        end do
    end do

end subroutine ffdev_err_pacpnl_error

! ==============================================================================
! subroutine ffdev_err_pacpnl_summary
! ==============================================================================

subroutine ffdev_err_pacpnl_summary(top,geo,printsum)

    use ffdev_targetset_dat
    use ffdev_geometry
    use ffdev_err_pacpnl_dat
    use ffdev_utils

    implicit none
    type(TOPOLOGY)  :: top
    type(GEOMETRY)  :: geo
    logical         :: printsum
    ! --------------------------------------------
    integer         :: k
    real(DEVDP)     :: totpnl, pnl, diff
    ! --------------------------------------------------------------------------

    if( .not. geo%sup_chrg_loaded) return

    if( printsum .eqv. .false. ) then
        printsum = .true.
        return
    end if

    write(DEV_OUT,*)
    write(DEV_OUT,5)
    write(DEV_OUT,10)
    write(DEV_OUT,20)

    totpnl = 0.0d0

    do k=1,top%natoms
        diff = geo%sup_chrg(k) - top%atoms(k)%charge
        pnl = 0.0

        select case(PACPnlFce)
            case(CHRG_PNL_QUADRATIC)
                pnl = PACPnlErrorWeight*diff**2
                totpnl = totpnl + pnl
            case default
                call ffdev_utils_exit(DEV_ERR,1,'Not implemented in ffdev_err_pacpnl_error!')
        end select

        write(DEV_OUT,30) k, top%atoms(k)%name, top%atom_types(top%atoms(k)%typeid)%name, top%atoms(k)%symmclass, &
                top%atoms(k)%charge, geo%sup_chrg(k), diff, pnl
    end do

    write(DEV_OUT,20)
    write(DEV_OUT,40) totpnl
    write(DEV_OUT,45) totpnl*PACPnlErrorWeight

 5 format('# Partial Atomic Charge Penalties')
10 format('# Atom  Name  Type SymmClass  TopCharge  GeoCharge Difference    Penalty')
20 format('# ---- ------ ---- --------- ---------- ---------- ---------- ----------')
30 format(I6,1X,A6,1X,A4,1X,I9,1X,F10.4,1X,F10.4,1X,F10.4,1X,F10.4)
40 format('# Final penalty               =                               ',F10.3)
45 format('# Final penalty (all weights) =                               ',F10.3)

end subroutine ffdev_err_pacpnl_summary

! ------------------------------------------------------------------------------

end module ffdev_err_pacpnl


