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
    SAPTRepErrorWeight     = 1.0
    SAPTDispErrorWeight    = 1.0

    SAPTErrorIndToEle      = .false.
    SAPTErrorIndToRep      = .true.
    SAPTErrorPenToRep      = .false.

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
    real(DEVDP)         :: err,serrele,serrrep,serrdis,trg_sapt_rep,pen_guess,totw
    real(DEVDP)         :: sapt_ele,trg_sapt_ele
    ! --------------------------------------------------------------------------

    error%sapt_ele = 0.0d0
    error%sapt_rep = 0.0d0
    error%sapt_dis = 0.0d0

    serrele = 0.0
    serrrep = 0.0
    serrdis = 0.0
    totw = 0.0

    do i=1,nsets
        if( .not. ( (sets(i)%nrefs .ge. 1) .or. (sets(i)%top%probe_size .gt. 0) ) ) cycle

        do j=1,sets(i)%ngeos
            ! ------------------------------------------------------------------
            if( .not. sets(i)%geo(j)%trg_sapt_loaded ) cycle

        ! electrostatics
            if( (.not. SAPTErrorPenToRep) .and. pen_enabled ) then
                sapt_ele = sets(i)%geo(j)%sapt_ele + sets(i)%geo(j)%sapt_pen
                trg_sapt_ele = sets(i)%geo(j)%trg_sapt_ele
                if( SAPTErrorIndToEle ) then
                    trg_sapt_ele = trg_sapt_ele + sets(i)%geo(j)%trg_sapt_ind
                end if
                err = sapt_ele - trg_sapt_ele
                serrele = serrele + sets(i)%geo(j)%weight * err**2
            end if

        ! repulsion
            trg_sapt_rep = sets(i)%geo(j)%trg_sapt_exc
            if( SAPTErrorIndToRep ) then
                trg_sapt_rep = trg_sapt_rep + sets(i)%geo(j)%trg_sapt_ind
            end if
            if( SAPTErrorPenToRep ) then
                ! pen_guess is penetration energy guess
                pen_guess = sets(i)%geo(j)%trg_sapt_ele - sets(i)%geo(j)%sapt_ele

                trg_sapt_rep = trg_sapt_rep + pen_guess
            end if
            err = sets(i)%geo(j)%sapt_rep - trg_sapt_rep
            serrrep = serrrep + sets(i)%geo(j)%weight * err**2

        ! dispersion
            err = sets(i)%geo(j)%sapt_dis - sets(i)%geo(j)%trg_sapt_dis
            serrdis = serrdis + sets(i)%geo(j)%weight * err**2

            totw = totw + sets(i)%geo(j)%weight
        end do
    end do

    if( totw .gt. 0 ) then
        error%sapt_ele = sqrt(serrele/totw)
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

    implicit none
    integer             :: i,j
    logical             :: printsum
    real(DEVDP)         :: serrele,serrrep,serrdisp,trg_sapt_rep,totw,sapt_ele
    real(DEVDP)         :: errele,errrep,errdisp,pen_guess,maxele,maxrep,maxdisp
    real(DEVDP)         :: trg_sapt_ele
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
    write(DEV_OUT,5)
    write(DEV_OUT,10)
    write(DEV_OUT,20)

    serrele = 0.0
    serrrep = 0.0
    serrdisp = 0.0

    maxele = 0.0d0
    maxrep = 0.0d0
    maxdisp = 0.0d0

    totw = 0.0d0

    do i=1,nsets
        printsum = .false.
        if( .not. ( (sets(i)%nrefs .ge. 1) .or. (sets(i)%top%probe_size .gt. 0) ) ) cycle

        do j=1,sets(i)%ngeos
            if( .not. sets(i)%geo(j)%trg_sapt_loaded ) cycle

            printsum = .true.

        ! electrostatics
            if( (.not. SAPTErrorPenToRep) .and. pen_enabled ) then
                sapt_ele = sets(i)%geo(j)%sapt_ele + sets(i)%geo(j)%sapt_pen
                trg_sapt_ele = sets(i)%geo(j)%trg_sapt_ele
                if( SAPTErrorIndToEle ) then
                    trg_sapt_ele = trg_sapt_ele + sets(i)%geo(j)%trg_sapt_ind
                end if
                errele = sapt_ele - trg_sapt_ele
                serrele = serrele + sets(i)%geo(j)%weight * errele**2
                if( abs(errele) .gt. abs(maxele) ) then
                    maxele = errele
                end if
            end if

            pen_guess = 0.0d0

        ! repulsion
            trg_sapt_rep = sets(i)%geo(j)%trg_sapt_exc
            if( SAPTErrorIndToRep ) then
                trg_sapt_rep = trg_sapt_rep + sets(i)%geo(j)%trg_sapt_ind
            end if
            if( SAPTErrorPenToRep ) then
                ! pen_guess is penetration energy guess
                pen_guess = sets(i)%geo(j)%trg_sapt_ele - sets(i)%geo(j)%sapt_ele
                trg_sapt_rep = trg_sapt_rep + pen_guess
            end if
            errrep = sets(i)%geo(j)%sapt_rep - trg_sapt_rep
            if( abs(errrep) .gt. abs(maxrep) ) then
                maxrep = errrep
            end if
            serrrep = serrrep + sets(i)%geo(j)%weight * errrep**2

        ! dispersion
            errdisp = sets(i)%geo(j)%sapt_dis - sets(i)%geo(j)%trg_sapt_dis
            if( abs(errdisp) .gt. abs(maxdisp) ) then
                maxdisp = errdisp
            end if
            serrdisp = serrdisp + sets(i)%geo(j)%weight * errdisp**2

            totw = totw + sets(i)%geo(j)%weight

            write(DEV_OUT,30) i, j, sets(i)%geo(j)%weight, &
                    sets(i)%geo(j)%trg_sapt_ele, sets(i)%geo(j)%sapt_ele, sets(i)%geo(j)%sapt_pen, errele,  &
                    pen_guess, sets(i)%geo(j)%trg_sapt_ind, sets(i)%geo(j)%trg_sapt_exc, &
                    trg_sapt_rep, sets(i)%geo(j)%sapt_rep, errrep, &
                    sets(i)%geo(j)%trg_sapt_dis, sets(i)%geo(j)%sapt_dis, errdisp
        end do
        if( printsum ) write(DEV_OUT,20)
    end do

    errele  = 0.0
    errrep  = 0.0
    errdisp = 0.0
    if( totw .gt. 0 ) then
        errele  = sqrt(serrele/totw)
        errrep  = sqrt(serrrep/totw)
        errdisp = sqrt(serrdisp/totw)
    end if

    write(DEV_OUT,35)  maxele, maxrep, maxdisp
    write(DEV_OUT,40)  errele, errrep, errdisp
    write(DEV_OUT,45)  SAPTEleErrorWeight*errele, SAPTRepErrorWeight*errrep, SAPTDispErrorWeight*errdisp

 5 format('# SAPT errors')
10 format('# SET GeoID Weight   ELE(TGR)    ELE(MM)    PEN(MM)   Err(ELE) PEN(GUESS)   IND(TGR)   EXC(TGR)   REP(TRG)' &
          '    REP(MM)   Err(REP)  DISP(TGR)   DISP(MM)  Err(DISP)')
20 format('# --- ----- ------ ---------- ---------- ---------- ---------- ---------- ---------- ---------- ----------' &
          ' ---------- ---------- ---------- ---------- ----------')
30 format(I5,1X,I5,1X,F6.3,1X,F10.3,1X,F10.3,1X,F10.3,1X,F10.3,1X,F10.3,1X,F10.3,1X,F10.3,1X,F10.3, &
         1X,F10.3,1X,F10.3,1X,F10.3,1X,F10.3,1X,F10.3)

35 format('# Maximum signed error (MSE)           = ',11X,F10.3,56X,F10.3,23X,F10.3)
40 format('# Root mean square error (RMSE)        = ',11X,F10.3,56X,F10.3,23X,F10.3)
45 format('# Final RMSE (all weights)             = ',11X,F10.3,56X,F10.3,23X,F10.3)

end subroutine ffdev_err_sapt_summary

! ------------------------------------------------------------------------------

end module ffdev_err_sapt


