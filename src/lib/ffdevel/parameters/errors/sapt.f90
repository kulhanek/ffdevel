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

    SAPTRepErrorWeight     = 1.0
    SAPTDispErrorWeight    = 1.0

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
    integer             :: i,j,nene
    real(DEVDP)         :: err,serrrep,serrdis,trg_sapt_rep,pen_guess
    ! --------------------------------------------------------------------------

    error%sapt_rep = 0.0d0
    error%sapt_dis = 0.0d0

    serrrep = 0.0
    serrdis = 0.0
    nene = 0

    do i=1,nsets
        if( sets(i)%isref ) cycle ! ignore reference sets

        do j=1,sets(i)%ngeos
            ! ------------------------------------------------------------------
            if( .not. (sets(i)%geo(j)%trg_sapt_loaded .and. sets(i)%nrefs .gt. 1) ) cycle

            nene = nene + 1

        ! electrostatics
            ! pen_guess is penetration energy guess
            pen_guess = sets(i)%geo(j)%trg_sapt_ele - sets(i)%geo(j)%sapt_ele

        ! repulsion
            trg_sapt_rep = sets(i)%geo(j)%trg_sapt_exc
            if( SAPTErrorIndToRep ) then
                trg_sapt_rep = trg_sapt_rep + sets(i)%geo(j)%trg_sapt_ind
            end if
            if( SAPTErrorPenToRep ) then
                trg_sapt_rep = trg_sapt_rep + pen_guess
            end if
            err = trg_sapt_rep - sets(i)%geo(j)%sapt_rep
            serrrep = serrrep + sets(i)%geo(j)%weight * err**2

        ! dispersion
            err = sets(i)%geo(j)%trg_sapt_dis - sets(i)%geo(j)%sapt_dis
            serrdis = serrdis + sets(i)%geo(j)%weight * err**2
        end do
    end do

    if( nene .gt. 0 ) then
        error%sapt_rep = sqrt(serrrep/real(nene))
        error%sapt_dis = sqrt(serrdis/real(nene))
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
    integer             :: i,j,nene
    logical             :: printsum
    real(DEVDP)         :: serrrep,serrdisp,trg_sapt_rep
    real(DEVDP)         :: errrep,errdisp,pen_guess
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

    serrrep = 0.0
    serrdisp = 0.0
    nene = 0

    do i=1,nsets
        printsum = .false.
        do j=1,sets(i)%ngeos
            if( .not. (sets(i)%geo(j)%trg_sapt_loaded .and. sets(i)%nrefs .gt. 1) ) cycle

            printsum = .true.
            nene = nene + 1

        ! electrostatics
            ! pen_guess is penetration energy guess
            pen_guess = sets(i)%geo(j)%trg_sapt_ele - sets(i)%geo(j)%sapt_ele

        ! repulsion
            trg_sapt_rep = sets(i)%geo(j)%trg_sapt_exc
            if( SAPTErrorIndToRep ) then
                trg_sapt_rep = trg_sapt_rep + sets(i)%geo(j)%trg_sapt_ind
            end if
            if( SAPTErrorPenToRep ) then
                trg_sapt_rep = trg_sapt_rep + pen_guess
            end if
            errrep = trg_sapt_rep - sets(i)%geo(j)%sapt_rep
            serrrep = serrrep + sets(i)%geo(j)%weight * errrep**2

        ! dispersion
            errdisp = sets(i)%geo(j)%trg_sapt_dis - sets(i)%geo(j)%sapt_dis
            serrdisp = serrdisp + sets(i)%geo(j)%weight * errdisp**2

            write(DEV_OUT,30) i, j, sets(i)%geo(j)%weight, &
                              sets(i)%geo(j)%trg_sapt_ele, sets(i)%geo(j)%sapt_ele, pen_guess, &
                              sets(i)%geo(j)%trg_sapt_ind, sets(i)%geo(j)%trg_sapt_exc, &
                              trg_sapt_rep, sets(i)%geo(j)%sapt_rep, errrep, &
                              sets(i)%geo(j)%trg_sapt_dis, sets(i)%geo(j)%sapt_dis, errdisp
        end do
        if( printsum ) write(DEV_OUT,20)
    end do

    errrep  = 0.0
    errdisp = 0.0
    if( nene .gt. 0 ) then
        errrep  = sqrt(serrrep/real(nene))
        errdisp = sqrt(serrdisp/real(nene))
    end if

    write(DEV_OUT,40)  errrep, errdisp
    write(DEV_OUT,45)  SAPTRepErrorWeight*errrep, SAPTDispErrorWeight*errdisp

 5 format('# SAPT errors')
10 format('# SET GeoID Weight   ELE(TGR)    ELE(MM) PEN(GUESS)   IND(TGR)   EXC(TGR)   REP(TRG)' &
          '    REP(MM) abs E(Err)  DISP(TGR)   DISP(MM) abs E(Err)')
20 format('# --- ----- ------ ---------- ---------- ---------- ---------- ---------- ----------' &
          ' ---------- ---------- ---------- ---------- ----------')
30 format(I5,1X,I5,1X,F6.3,1X,F10.3,1X,F10.3,1X,F10.3,1X,F10.3,1X,F10.3,1X,F10.3,&
         1X,F10.3,1X,F10.3,1X,F10.3,1X,F10.3,1X,F10.3)

40 format('# Final error (weighted per geometry)  = ',55X,F10.3,23X,F10.3)
45 format('# Final error (all weights)            = ',55X,F10.3,23X,F10.3)

end subroutine ffdev_err_sapt_summary

! ------------------------------------------------------------------------------

end module ffdev_err_sapt

