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

module ffdev_err_sapt

use ffdev_constants
use ffdev_variables

contains

! ==============================================================================
! subroutine ffdev_err_sapt_init
! ==============================================================================

subroutine ffdev_err_sapt_init

    use ffdev_err_sapt_dat
    use ffdev_errors_dat

    implicit none
    ! --------------------------------------------------------------------------

    EnableSAPTError        = .false.
    PrintSAPTErrorSummary  = .false.

    SAPTEleErrorWeight     = 1.0
    SAPTIndErrorWeight     = 1.0
    SAPTRepErrorWeight     = 1.0
    SAPTDisErrorWeight    = 1.0

    SAPTErrorIndToRep      = .true.
    SAPTErrorPenToRep      = .true.

end subroutine ffdev_err_sapt_init

! ==============================================================================
! subroutine ffdev_err_sapt_error
! ==============================================================================

subroutine ffdev_err_sapt_error(error)

    use ffdev_targetset
    use ffdev_targetset_dat
    use ffdev_utils
    use ffdev_geometry
    use ffdev_errors_dat
    use ffdev_err_sapt_dat

    implicit none
    type(FFERROR_TYPE)  :: error
    ! --------------------------------------------
    integer             :: i,j
    real(DEVDP)         :: err,serrele,serrrep,serrdis,pen_guess,totw
    real(DEVDP)         :: sapt_ele,trg_sapt_exc,serrind
    ! --------------------------------------------------------------------------

    error%sapt_ele = 0.0d0
    error%sapt_ind = 0.0d0
    error%sapt_rep = 0.0d0
    error%sapt_dis = 0.0d0

    serrele = 0.0
    serrind = 0.0
    serrrep = 0.0
    serrdis = 0.0
    totw = 0.0

    do i=1,nsets
        if( .not. ( (sets(i)%nrefs .ge. 1) .or. (sets(i)%top%probe_size .gt. 0) ) ) cycle

        do j=1,sets(i)%ngeos
            ! ------------------------------------------------------------------
            if( .not. sets(i)%geo(j)%trg_sapt_loaded ) cycle

        ! electrostatics
            if( pen_enabled ) then
                sapt_ele = sets(i)%geo(j)%sapt_ele + sets(i)%geo(j)%sapt_pen
                err = sapt_ele - sets(i)%geo(j)%trg_sapt_ele
                serrele = serrele + sets(i)%geo(j)%weight * err**2
            end if

        ! induction
            if( ind_enabled ) then
                err = sets(i)%geo(j)%sapt_ind - sets(i)%geo(j)%trg_sapt_ind
                serrind = serrind + sets(i)%geo(j)%weight * err**2
            end if

        ! repulsion
            trg_sapt_exc = sets(i)%geo(j)%trg_sapt_exc
            if( SAPTErrorIndToRep .and. (.not. ind_enabled)  ) then
                trg_sapt_exc = trg_sapt_exc + sets(i)%geo(j)%trg_sapt_ind
            end if
            if( SAPTErrorPenToRep .and. (.not. pen_enabled)  ) then
                ! pen_guess is penetration energy guess
                pen_guess = sets(i)%geo(j)%trg_sapt_ele - sets(i)%geo(j)%sapt_ele
                trg_sapt_exc = trg_sapt_exc + pen_guess
            end if
            err = sets(i)%geo(j)%sapt_rep - trg_sapt_exc
            serrrep = serrrep + sets(i)%geo(j)%weight * err**2

        ! dispersion
            err = sets(i)%geo(j)%sapt_dis - sets(i)%geo(j)%trg_sapt_dis
            serrdis = serrdis + sets(i)%geo(j)%weight * err**2

            totw = totw + sets(i)%geo(j)%weight
        end do
    end do

    if( totw .gt. 0 ) then
        error%sapt_ele = sqrt(serrele/totw)
        error%sapt_ind = sqrt(serrind/totw)
        error%sapt_rep = sqrt(serrrep/totw)
        error%sapt_dis = sqrt(serrdis/totw)
    end if

end subroutine ffdev_err_sapt_error

! ==============================================================================
! subroutine ffdev_err_sapt_summary
! ==============================================================================

subroutine ffdev_err_sapt_summary

    use ffdev_targetset_dat
    use ffdev_geometry
    use ffdev_err_sapt_dat
    use prmfile

    implicit none
    integer             :: i,j
    logical             :: printsum
    real(DEVDP)         :: serrele,serrind,serrrep,serrdis,totw
    real(DEVDP)         :: errele,errind,errrep,errdis,pen_guess
    real(DEVDP)         :: maxele,maxind,maxrep,maxdis,sapt_ele,trg_sapt_exc
    ! --------------------------------------------------------------------------

    printsum = .false.
    do i=1,nsets
        do j=1,sets(i)%ngeos
            if( sets(i)%geo(j)%trg_sapt_loaded .and. (sets(i)%isref .eqv. .false.) ) then
                printsum = .true.
            end if
        end do
    end do
    if( .not. printsum ) return

    write(DEV_OUT,*)
    write(DEV_OUT,5) trim(prmfile_onoff(SAPTErrorPenToRep .and. (.not. pen_enabled))), &
                     trim(prmfile_onoff(SAPTErrorIndToRep .and. (.not. ind_enabled)))
    write(DEV_OUT, 10,ADVANCE='NO')
    write(DEV_OUT,110,ADVANCE='NO')
    write(DEV_OUT,210,ADVANCE='NO')
    write(DEV_OUT,310,ADVANCE='NO')
    write(DEV_OUT,410)

    write(DEV_OUT, 20,ADVANCE='NO')
    write(DEV_OUT,120,ADVANCE='NO')
    write(DEV_OUT,220,ADVANCE='NO')
    write(DEV_OUT,320,ADVANCE='NO')
    write(DEV_OUT,420)

    serrele = 0.0
    serrind = 0.0
    serrrep = 0.0
    serrdis = 0.0

    maxele = 0.0d0
    maxind = 0.0d0
    maxrep = 0.0d0
    maxdis = 0.0d0

    totw = 0.0d0

    do i=1,nsets
        printsum = .false.
        if( .not. ( (sets(i)%nrefs .ge. 1) .or. (sets(i)%top%probe_size .gt. 0) ) ) cycle

        do j=1,sets(i)%ngeos
            if( .not. sets(i)%geo(j)%trg_sapt_loaded ) cycle

            printsum = .true.

        ! electrostatics
            errele = 0.0d0
            sapt_ele = 0.0d0
            if( pen_enabled ) then
                sapt_ele = sets(i)%geo(j)%sapt_ele + sets(i)%geo(j)%sapt_pen
                errele = sapt_ele - sets(i)%geo(j)%trg_sapt_ele
                serrele = serrele + sets(i)%geo(j)%weight * errele**2
                if( abs(errele) .gt. abs(maxele) ) then
                    maxele = errele
                end if
            end if

        ! induction
            errind = 0.0d0
            if( ind_enabled ) then
                errind = sets(i)%geo(j)%sapt_ind - sets(i)%geo(j)%trg_sapt_ind
                serrind = serrind + sets(i)%geo(j)%weight * errind**2
                if( abs(errind) .gt. abs(maxind) ) then
                    maxind = errind
                end if
            end if

        ! repulsion
            pen_guess = 0.0d0
            trg_sapt_exc = sets(i)%geo(j)%trg_sapt_exc
            if( SAPTErrorIndToRep .and. (.not. ind_enabled) ) then
                trg_sapt_exc = trg_sapt_exc + sets(i)%geo(j)%trg_sapt_ind
            end if
            if( SAPTErrorPenToRep .and. (.not. pen_enabled) ) then
                ! pen_guess is penetration energy guess
                pen_guess = sets(i)%geo(j)%trg_sapt_ele - sets(i)%geo(j)%sapt_ele
                trg_sapt_exc = trg_sapt_exc + pen_guess
            end if
            errrep = sets(i)%geo(j)%sapt_rep - trg_sapt_exc
            if( abs(errrep) .gt. abs(maxrep) ) then
                maxrep = errrep
            end if
            serrrep = serrrep + sets(i)%geo(j)%weight * errrep**2

        ! dispersion
            errdis = sets(i)%geo(j)%sapt_dis - sets(i)%geo(j)%trg_sapt_dis
            if( abs(errdis) .gt. abs(maxdis) ) then
                maxdis = errdis
            end if
            serrdis = serrdis + sets(i)%geo(j)%weight * errdis**2

            totw = totw + sets(i)%geo(j)%weight

            write(DEV_OUT, 30,ADVANCE='NO') i, j, sets(i)%geo(j)%weight
            write(DEV_OUT,130,ADVANCE='NO') sets(i)%geo(j)%sapt_ele, sets(i)%geo(j)%sapt_pen, &
                                            sapt_ele, sets(i)%geo(j)%trg_sapt_ele, errele
            write(DEV_OUT,230,ADVANCE='NO') sets(i)%geo(j)%sapt_ind, sets(i)%geo(j)%trg_sapt_ind, errind
            write(DEV_OUT,330,ADVANCE='NO') pen_guess, sets(i)%geo(j)%sapt_rep, &
                                            sets(i)%geo(j)%trg_sapt_exc, trg_sapt_exc, errrep
            write(DEV_OUT,430)              sets(i)%geo(j)%sapt_dis, sets(i)%geo(j)%trg_sapt_dis, errdis
        end do
        if( printsum ) then
            write(DEV_OUT, 20,ADVANCE='NO')
            write(DEV_OUT,120,ADVANCE='NO')
            write(DEV_OUT,220,ADVANCE='NO')
            write(DEV_OUT,320,ADVANCE='NO')
            write(DEV_OUT,420)
        end if
    end do

    errele = 0.0
    errind = 0.0
    errrep = 0.0
    errdis = 0.0
    if( totw .gt. 0 ) then
        errele = sqrt(serrele/totw)
        errind = sqrt(serrind/totw)
        errrep = sqrt(serrrep/totw)
        errdis = sqrt(serrdis/totw)
    end if

    write(DEV_OUT,35)  maxele, maxind, maxrep, maxdis
    write(DEV_OUT,40)  errele, errind, errrep, errdis
    write(DEV_OUT,45)  SAPTEleErrorWeight*errele, SAPTIndErrorWeight*errind,  &
                       SAPTRepErrorWeight*errrep, SAPTDisErrorWeight*errdis

  5 format('# SAPT errors with effective setup: penasrep = ',A,'; indasrep = ',A)
 10 format('# SET GeoID Weight | ')
 20 format('# --- ----- ------ | ')
 30 format(I5,1X,I5,1X,F6.3,3X)

110 format('   ELE(MM)    PEN(MM) TOTELE(MM)   ELE(TRG)   Err(ELE) | ')
120 format('---------- ---------- ---------- ---------- ---------- | ')
130 format(F10.3,1X,F10.3,1X,F10.3,1X,F10.3,1X,F10.3,3X)

210 format('   IND(MM)   IND(TRG)   Err(IND) | ')
220 format('---------- ---------- ---------- | ')
230 format(F10.3,1X,F10.3,1X,F10.3,3X)

310 format('PEN(GUESS)    REP(MM)   EXC(TRG)   REP(TRG)   Err(REP) | ')
320 format('---------- ---------- ---------- ---------- ---------- | ')
330 format(F10.3,1X,F10.3,1X,F10.3,1X,F10.3,1X,F10.3,3X)

410 format('   DIS(MM)   DIS(TRG)   Err(DIS)')
420 format('---------- ---------- ----------')
430 format(F10.3,1X,F10.3,1X,F10.3)

!          '# SET GeoID Weight |    ELE(MM)    PEN(MM) TOTELE(MM)   ELE(TRG)   Err(ELE) | '
!          '                     ---------- ---------- ---------- ---------- ---------- | '
 35 format('# Maximum signed error (MSE)             = ',22X,F10.3,25X,F10.3,47X,F10.3,25X,F10.3)
 40 format('# Root mean square error (RMSE)          = ',22X,F10.3,25X,F10.3,47X,F10.3,25X,F10.3)
 45 format('# Final RMSE (all weights)               = ',22X,F10.3,25X,F10.3,47X,F10.3,25X,F10.3)

end subroutine ffdev_err_sapt_summary

! ------------------------------------------------------------------------------

end module ffdev_err_sapt


