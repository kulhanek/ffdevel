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
    PACPnlErrorTempFactor    = 0.0

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
    real(DEVDP)         :: diff,totw,tote,w,molerr
    ! --------------------------------------------------------------------------

    error%pacpnl = 0.0d0

    tote = 0.0d0
    totw = 0.0d0

    do i=1,nsets
        do j=1,sets(i)%ngeos
            if( .not. sets(i)%geo(j)%sup_chrg_loaded ) cycle
            w = 1.0d0
            if( (PACPnlErrorTempFactor .gt. 0) .and. sets(i)%geo(j)%trg_ene_loaded ) then
                w = exp(-sets(i)%geo(j)%trg_energy/(PACPnlErrorTempFactor*DEV_Rgas))
            end if
            molerr = 0.0
            do k=1,sets(i)%top%natoms
                select case(PACPnlFce)
                    case(CHRG_PNL_QUADRATIC)
                        diff = sets(i)%top%atoms(k)%charge - sets(i)%geo(j)%sup_chrg(k)
                        molerr = molerr + diff**2
                    case default
                        call ffdev_utils_exit(DEV_ERR,1,'Not implemented in ffdev_err_pacpnl_error!')
                end select
            end do
            if( sets(i)%top%natoms .gt. 0 ) then
                tote = tote + w*sqrt(molerr/real(sets(i)%top%natoms,DEVDP))
                totw = totw + w
            end if
        end do
    end do

    if( totw .gt. 0.0d0 ) then
        error%pacpnl = tote/totw
    end if

end subroutine ffdev_err_pacpnl_error

! ==============================================================================
! subroutine ffdev_err_pacpnl_summary
! ==============================================================================

subroutine ffdev_err_pacpnl_summary

    use ffdev_targetset_dat
    use ffdev_geometry
    use ffdev_err_pacpnl_dat
    use ffdev_utils

    implicit none
    integer         :: i,j,k
    logical         :: printsum
    real(DEVDP)     :: pnl, diff, w, molerr, tote, totw, maxv, minv
    ! --------------------------------------------------------------------------

    printsum = .false.
    do i=1,nsets
        do j=1,sets(i)%ngeos
            if( .not. sets(i)%geo(j)%sup_chrg_loaded ) cycle
            printsum = .true.
            exit
        end do
    end do

    write(DEV_OUT,*)
    write(DEV_OUT,5)
    write(DEV_OUT,110)
    write(DEV_OUT,120)

    tote = 0.0d0
    totw = 0.0d0

    do i=1,nsets
        printsum = .false.
        do j=1,sets(i)%ngeos
            if( .not. sets(i)%geo(j)%sup_chrg_loaded ) cycle
            w = 1.0d0
            if( (PACPnlErrorTempFactor .gt. 0) .and. sets(i)%geo(j)%trg_ene_loaded ) then
                w = exp(-sets(i)%geo(j)%trg_energy/(PACPnlErrorTempFactor*DEV_Rgas))
            end if
            molerr = 0.0
            do k=1,sets(i)%top%natoms
                select case(PACPnlFce)
                    case(CHRG_PNL_QUADRATIC)
                        diff = sets(i)%top%atoms(k)%charge - sets(i)%geo(j)%sup_chrg(k)
                        molerr = molerr + diff**2
                    case default
                        call ffdev_utils_exit(DEV_ERR,1,'Not implemented in ffdev_err_pacpnl_error!')
                end select
            end do

            pnl = 0.0d0
            if( sets(i)%top%natoms .gt. 0 ) then
                pnl = sqrt(molerr/real(sets(i)%top%natoms,DEVDP))
                tote = tote + w*pnl
                totw = totw + w
                printsum = .true.
                write(DEV_OUT,130) i,j,w,pnl
            end if

        end do
        if( printsum ) write(DEV_OUT,120)
    end do

    if( totw .gt. 0.0d0 ) then
        tote = tote/totw
    end if

    write(DEV_OUT,140) tote
    write(DEV_OUT,150) tote*PACPnlErrorWeight

! ----------------------------

    write(DEV_OUT,*)
    write(DEV_OUT,6)

    do i=1,nsets
        do j=1,sets(i)%ngeos
            if( .not. sets(i)%geo(j)%sup_chrg_loaded ) cycle

            write(DEV_OUT,*)
            write(DEV_OUT,7) i,j

            write(DEV_OUT,10)
            write(DEV_OUT,20)

            maxv = 0.0d0
            minv = 0.0d0

            do k=1,sets(i)%top%natoms
                diff = sets(i)%geo(j)%sup_chrg(k) - sets(i)%top%atoms(k)%charge

                write(DEV_OUT,30) k, sets(i)%top%atoms(k)%name, sets(i)%top%atom_types(sets(i)%top%atoms(k)%typeid)%name, &
                        sets(i)%top%atoms(k)%symmclass, &
                        sets(i)%top%atoms(k)%charge, sets(i)%geo(j)%sup_chrg(k), diff

                if( k .eq. 1 ) then
                    maxv = diff
                    minv = diff
                end if
                if( maxv .lt. diff ) then
                    maxv = diff
                end if
                if( minv .gt. diff ) then
                    minv = diff
                end if
            end do

            write(DEV_OUT,20)
            write(DEV_OUT,60) minv
            write(DEV_OUT,70) maxv

        end do
    end do

  5 format('# Partial Atomic Charge Penalties')
  6 format('# Partial Atomic Charge Summaries Per molecules')

110 format('# SET  GeoID Weight    Penalty')
120 format('# --- ------ ------ ----------')
130 format(I5,1X,I6,1X,F6.4,1X,F10.4)
140 format('# Final error (weighted per geometry) =  ',F10.3)
150 format('# Final error (all weights)           =  ',F10.3)

  7 format('== [SET ',I5.5,']/[GEO ',I6.6,'] ====================================================')

 10 format('# Atom  Name  Type SymmClass  TopCharge  GeoCharge Difference')
 20 format('# ---- ------ ---- --------- ---------- ---------- ----------')
 30 format(I6,1X,A6,1X,A4,1X,I9,1X,F10.4,1X,F10.4,1X,F10.4,1X,F10.4,1X,F10.4)
 60 format('# Min difference       =                           ',F10.4)
 70 format('# Max difference       =                           ',F10.4)

end subroutine ffdev_err_pacpnl_summary

! ------------------------------------------------------------------------------

end module ffdev_err_pacpnl


