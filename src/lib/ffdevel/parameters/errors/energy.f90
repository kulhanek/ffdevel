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

module ffdev_err_energy

use ffdev_constants
use ffdev_variables

contains

! ==============================================================================
! subroutine ffdev_err_energy_init
! ==============================================================================

subroutine ffdev_err_energy_init

    use ffdev_err_energy_dat
    use ffdev_errors_dat

    implicit none
    ! --------------------------------------------------------------------------

    EnableEnergyError       = .false.
    PrintEnergyErrorSummary = .false.
    EnergyErrorWeight       = 1.0
    EnergyErrorMode         = EE_ABS
    EnableMaxFilter         = .false.
    MaxTargetEnergy         = 0.0
    EnableMinFilter         = .false.
    MinTargetEnergy         = 0.0

end subroutine ffdev_err_energy_init

! ==============================================================================
! subroutine ffdev_err_energy_error
! ==============================================================================

subroutine ffdev_err_energy_error(error)

    use ffdev_targetset
    use ffdev_targetset_dat
    use ffdev_utils
    use ffdev_geometry
    use ffdev_errors_dat
    use ffdev_err_energy_dat

    implicit none
    type(FFERROR_TYPE)  :: error
    ! --------------------------------------------
    integer             :: i,j,nene
    real(DEVDP)         :: err,seterrene,totw
    ! --------------------------------------------------------------------------

    error%energy = 0.0d0

    seterrene = 0.0
    nene = 0
    totw = 0

    do i=1,nsets
        ! use only sets, which can provide reliable energy
        if( .not. ( (sets(i)%nrefs .ge. 1) .or. (sets(i)%top%probe_size .gt. 0) ) ) cycle

        do j=1,sets(i)%ngeos
            ! ------------------------------------------------------------------
            if( .not. sets(i)%geo(j)%trg_ene_loaded ) cycle

            ! filters
            if( EnableMaxFilter ) then
                if( sets(i)%geo(j)%trg_energy .gt. MaxTargetEnergy ) cycle
            end if
            if( EnableMinFilter ) then
                if( sets(i)%geo(j)%trg_energy .lt. MinTargetEnergy ) cycle
            end if

            select case(EnergyErrorMode)
                case(EE_ABS)
                    nene = nene + 1
                    err = sets(i)%geo(j)%total_ene - sets(i)%geo(j)%trg_energy
                    seterrene = seterrene + sets(i)%geo(j)%weight * err**2
                    totw = totw + sets(i)%geo(j)%weight
                case(EE_REL)
                    if( sets(i)%geo(j)%trg_energy .gt. 0 ) then
                        nene = nene + 1
                        err = (sets(i)%geo(j)%total_ene - sets(i)%geo(j)%trg_energy) / &
                              sets(i)%geo(j)%trg_energy
                        seterrene = seterrene + sets(i)%geo(j)%weight * err**2
                        totw = totw + sets(i)%geo(j)%weight
                    end if
                case(EE_LOG)
                    if( (sets(i)%geo(j)%total_ene .gt. 0) .and. &
                        (sets(i)%geo(j)%trg_energy .gt. 0) ) then
                        nene = nene + 1
                        err = log(sets(i)%geo(j)%total_ene) - log(sets(i)%geo(j)%trg_energy)
                        seterrene = seterrene + sets(i)%geo(j)%weight * err**2
                        totw = totw + sets(i)%geo(j)%weight
                    end if
            end select

        end do
    end do

    if( totw .gt. 0 ) then
        error%energy = sqrt(seterrene/totw)
    end if

end subroutine ffdev_err_energy_error

! ==============================================================================
! subroutine ffdev_err_energy_summary
! ==============================================================================

subroutine ffdev_err_energy_summary

    use ffdev_targetset_dat
    use ffdev_geometry
    use ffdev_err_energy_dat

    implicit none
    real(DEVDP)         :: aerr,aserr
    real(DEVDP)         :: rerr,rserr
    real(DEVDP)         :: lerr,lserr,maxerr,atotw,rtotw,ltotw
    integer             :: i,j
    logical             :: printsum
    ! --------------------------------------------------------------------------

    printsum = .false.
    do i=1,nsets
        do j=1,sets(i)%ngeos
            if( sets(i)%geo(j)%trg_ene_loaded .and. (sets(i)%isref .eqv. .false.) ) then
                printsum = .true.
            end if
        end do
    end do
    if( .not. printsum ) return

    write(DEV_OUT,*)
    write(DEV_OUT,5)
    write(DEV_OUT,10,ADVANCE='NO')
    write(DEV_OUT,50)
    write(DEV_OUT,20,ADVANCE='NO')
    write(DEV_OUT,60)

    aserr = 0.0d0
    rserr = 0.0d0
    lserr = 0.0d0
    atotw = 0.0d0
    rtotw = 0.0d0
    ltotw = 0.0d0

    maxerr = 0.0d0

    do i=1,nsets
        ! use only sets, which can provide reliable energy
        if( .not. ( (sets(i)%nrefs .ge. 1) .or. (sets(i)%top%probe_size .gt. 0) ) ) cycle

        do j=1,sets(i)%ngeos
            printsum = .false.
            if( .not. sets(i)%geo(j)%trg_ene_loaded ) cycle

            printsum = .true.

            ! absolute
            aerr  = sets(i)%geo(j)%total_ene - sets(i)%geo(j)%trg_energy
            aserr = aserr + sets(i)%geo(j)%weight * aerr**2
            atotw  = atotw + sets(i)%geo(j)%weight

            if( abs(aerr) .gt. abs(maxerr) ) then
                maxerr = aerr
            end if

            ! relative
            rerr  = 0.0
            if( sets(i)%geo(j)%trg_energy .ne. 0 ) then
                rerr  = (sets(i)%geo(j)%total_ene - sets(i)%geo(j)%trg_energy)/sets(i)%geo(j)%trg_energy
                rserr = rserr + sets(i)%geo(j)%weight * rerr**2
                rtotw  = rtotw + sets(i)%geo(j)%weight
            end if
            ! log
            lerr  = 0.0
            if( (sets(i)%geo(j)%trg_energy .gt. 0) .and. (sets(i)%geo(j)%total_ene .gt. 0) ) then
                lerr  = log(sets(i)%geo(j)%total_ene) - log(sets(i)%geo(j)%trg_energy)
                lserr = lserr + sets(i)%geo(j)%weight * lerr**2
                ltotw  = ltotw + sets(i)%geo(j)%weight
            end if

            write(DEV_OUT,30,ADVANCE='NO') i, j, sets(i)%geo(j)%weight, &
                              sets(i)%geo(j)%total_ene, sets(i)%geo(j)%trg_energy, aerr, rerr*100.0d0, lerr

            write(DEV_OUT,70) sets(i)%geo(j)%bond_ene, sets(i)%geo(j)%angle_ene, sets(i)%geo(j)%dih_ene, &
                              sets(i)%geo(j)%impropr_ene, sets(i)%geo(j)%dih_ene + sets(i)%geo(j)%impropr_ene, &
                              sets(i)%geo(j)%ele_ene, sets(i)%geo(j)%pen_ene, sets(i)%geo(j)%ele14_ene, &
                              sets(i)%geo(j)%ele_ene + sets(i)%geo(j)%pen_ene + sets(i)%geo(j)%ele14_ene, &
                              sets(i)%geo(j)%ind_ene, &
                              sets(i)%geo(j)%rep_ene, sets(i)%geo(j)%rep14_ene, sets(i)%geo(j)%rep_ene + sets(i)%geo(j)%rep14_ene, &
                              sets(i)%geo(j)%dis_ene, sets(i)%geo(j)%dis14_ene, sets(i)%geo(j)%dis_ene + sets(i)%geo(j)%dis14_ene, &
                              sets(i)%geo(j)%bond_ene + sets(i)%geo(j)%angle_ene      &
                               + sets(i)%geo(j)%dih_ene + sets(i)%geo(j)%impropr_ene, &
                              sets(i)%geo(j)%ele_ene + sets(i)%geo(j)%pen_ene + sets(i)%geo(j)%ele14_ene  &
                               + sets(i)%geo(j)%ind_ene &
                               + sets(i)%geo(j)%rep_ene + sets(i)%geo(j)%rep14_ene  &
                               + sets(i)%geo(j)%dis_ene + sets(i)%geo(j)%dis14_ene

            if( Verbosity .ge. DEV_VERBOSITY_FULL ) then
                call ffdev_geometry_info_ene(sets(i)%geo(j))
                write(DEV_OUT,*)
            end if

        end do
        if( printsum ) then
            write(DEV_OUT,20,ADVANCE='NO')
            write(DEV_OUT,60)
        end if
    end do

    if( atotw .gt. 0 ) then
        aserr = sqrt(aserr / atotw)
    end if
    if( rtotw .gt. 0 ) then
        rserr = sqrt(rserr / rtotw)
    end if
    if( ltotw .gt. 0 ) then
        lserr = sqrt(lserr / ltotw)
    end if

    write(DEV_OUT,35)  maxerr
    write(DEV_OUT,40)  aserr, rserr*100.0d0, lserr
    write(DEV_OUT,45)  EnergyErrorWeight*aserr, EnergyErrorWeight*rserr*100.0d0, EnergyErrorWeight*lserr

 5 format('# Energy errors')
10 format('# SET GeoID Weight      E(MM)     E(TGR)     Err(E) relErr%(E)  logErr(E) | ')
20 format('# --- ----- ------ ---------- ---------- ---------- ---------- ---------- | ')
30 format(I5,1X,I5,1X,F6.3,1X,F10.3,1X,F10.3,1X,F10.3,1X,F10.3,1X,F10.3,3X)

35 format('# Maximum signed error (MSE)          =  ',F10.3)
40 format('# Root mean square error (RMSE)       =  ',F10.3,1X,F10.3,1X,F10.3)
45 format('# Final RMSE (all weights)            =  ',F10.3,1X,F10.3,1X,F10.3)

50 format('     Ebonds    Eangles      Etors      Eimps  Edih(t+i)        Eel       Epen      E14el    Etotele' &
          '       Eind       Erep     E14rep    Etotrep      Edisp    E14disp   Etotdisp        Ebn        Enb')
60 format(' ---------- ---------- ---------- ---------- ---------- ---------- ---------- ---------- ----------' &
          ' ---------- ---------- ---------- ---------- ---------- ---------- ---------- ---------- ----------')
70 format(1X,F10.3,1X,F10.3,1X,F10.3,1X,F10.3,1X,F10.3,1X,F10.3,1X,F10.3,1X,F10.3,1X,F10.3, &
          1X,F10.3,1X,F10.3,1X,F10.3,1X,F10.3,1X,F10.3,1X,F10.3,1X,F10.3,1X,F10.3,1X,F10.3)

end subroutine ffdev_err_energy_summary

! ------------------------------------------------------------------------------

end module ffdev_err_energy


