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

    SAPTErrorIndToRep      = .true.

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
    real(DEVDP)         :: err,serrele,serrrep,serrdisp,trg_sapt_rep,trg_sapt_ele
    ! --------------------------------------------------------------------------

    error%sapt_ele = 0.0d0
    error%sapt_rep = 0.0d0
    error%sapt_disp = 0.0d0

    serrele = 0.0
    serrrep = 0.0
    serrdisp = 0.0
    nene = 0

    do i=1,nsets
        if( sets(i)%isref ) cycle ! ignore reference sets

        do j=1,sets(i)%ngeos
            ! ------------------------------------------------------------------
            if( .not. (sets(i)%geo(j)%trg_sapt_loaded .and. sets(i)%nrefs .gt. 1) ) cycle

            nene = nene + 1

        ! electrostatics
            trg_sapt_ele = sets(i)%geo(j)%trg_sapt_ele
            ! FIXME
!            trg_sapt_ele = trg_sapt_ele + sets(i)%geo(j)%trg_sapt_ind
            err = sets(i)%geo(j)%sapt_ele - trg_sapt_ele
            serrele = serrele + sets(i)%geo(j)%weight * err**2

        ! repulsion
            trg_sapt_rep = sets(i)%geo(j)%trg_sapt_exch
            ! FIXME
!            if( SAPTErrorIndToRep ) then
                trg_sapt_rep = trg_sapt_rep + sets(i)%geo(j)%trg_sapt_ind
!            end if
            err = sets(i)%geo(j)%sapt_rep - trg_sapt_rep
            serrrep = serrrep + sets(i)%geo(j)%weight * err**2

        ! dispersion
            err = sets(i)%geo(j)%sapt_disp - sets(i)%geo(j)%trg_sapt_disp
            serrdisp = serrdisp + sets(i)%geo(j)%weight * err**2
        end do
    end do

    if( nene .gt. 0 ) then
        error%sapt_ele  = sqrt(serrele/real(nene))
        error%sapt_rep  = sqrt(serrrep/real(nene))
        error%sapt_disp = sqrt(serrdisp/real(nene))
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
    real(DEVDP)         :: serrele,serrrep,serrdisp,trg_sapt_rep
    real(DEVDP)         :: errele,errrep,errdisp
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
    nene = 0

    do i=1,nsets
        printsum = .false.
        do j=1,sets(i)%ngeos
            if( .not. (sets(i)%geo(j)%trg_sapt_loaded .and. sets(i)%nrefs .gt. 1) ) cycle

            printsum = .true.
            nene = nene + 1

        ! electrostatics
            errele = sets(i)%geo(j)%sapt_ele - sets(i)%geo(j)%trg_sapt_ele
            serrele = serrele + sets(i)%geo(j)%weight * errele**2

        ! repulsion
            trg_sapt_rep = sets(i)%geo(j)%trg_sapt_exch
            if( SAPTErrorIndToRep ) then
                trg_sapt_rep = trg_sapt_rep + sets(i)%geo(j)%trg_sapt_ind
            end if
            errrep = sets(i)%geo(j)%sapt_rep - trg_sapt_rep
            serrrep = serrrep + sets(i)%geo(j)%weight * errrep**2

        ! dispersion
            errdisp = sets(i)%geo(j)%sapt_disp - sets(i)%geo(j)%trg_sapt_disp
            serrdisp = serrdisp + sets(i)%geo(j)%weight * errdisp**2


            write(DEV_OUT,30) i, j, sets(i)%geo(j)%weight, &
                              sets(i)%geo(j)%trg_sapt_ele, sets(i)%geo(j)%sapt_ele, errele, &
                              trg_sapt_rep, sets(i)%geo(j)%sapt_rep, errrep, &
                              sets(i)%geo(j)%trg_sapt_disp, sets(i)%geo(j)%sapt_disp, errdisp
        end do
        if( printsum ) write(DEV_OUT,20)
    end do

    errele  = 0.0
    errrep  = 0.0
    errdisp = 0.0
    if( nene .gt. 0 ) then
        errele  = sqrt(serrele/real(nene))
        errrep  = sqrt(serrrep/real(nene))
        errdisp = sqrt(serrdisp/real(nene))
    end if

    write(DEV_OUT,40)  errele, errrep, errdisp
    write(DEV_OUT,45)  SAPTEleErrorWeight*errele, SAPTRepErrorWeight*errrep, SAPTDispErrorWeight*errdisp

 5 format('# SAPT errors')
10 format('# SET GeoID Weight   ELE(TGR)    ELE(MM) abs E(Err)   REP(TGR)    REP(MM) abs E(Err)  DISP(TGR)   DISP(MM) abs E(Err)')
20 format('# --- ----- ------ ---------- ---------- ---------- ---------- ---------- ---------- ---------- ---------- ----------')
30 format(I5,1X,I5,1X,F6.3,1X,F10.3,1X,F10.3,1X,F10.3,1X,F10.3,1X,F10.3,1X,F10.3,1X,F10.3,1X,F10.3,1X,F10.3)

40 format('# Final error (weighted per geometry)  = ',F10.3,23X,F10.3,23X,F10.3)
45 format('# Final error (all weights)            = ',F10.3,23X,F10.3,23X,F10.3)

end subroutine ffdev_err_sapt_summary

! ------------------------------------------------------------------------------

end module ffdev_err_sapt


