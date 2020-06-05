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

module ffdev_err_probe

use ffdev_constants
use ffdev_variables

contains

! ==============================================================================
! subroutine ffdev_err_probe_init
! ==============================================================================

subroutine ffdev_err_probe_init

    use ffdev_err_probe_dat
    use ffdev_errors_dat

    implicit none
    ! --------------------------------------------------------------------------

    EnableProbeError        = .false.
    PrintProbeErrorSummary  = .false.

    ProbeErrorWeight        = 1.0

end subroutine ffdev_err_probe_init

! ==============================================================================
! subroutine ffdev_err_probe_error
! ==============================================================================

subroutine ffdev_err_probe_error(error)

    use ffdev_targetset
    use ffdev_targetset_dat
    use ffdev_utils
    use ffdev_geometry
    use ffdev_errors_dat
    use ffdev_err_probe_dat

    implicit none
    type(FFERROR_TYPE)  :: error
    ! --------------------------------------------
    integer             :: i,j,nene
    real(DEVDP)         :: err,serene,emm
    ! --------------------------------------------------------------------------

    error%probe_ene = 0.0d0

    serene = 0.0
    nene = 0

    do i=1,nsets
        if( sets(i)%top%probe_size .eq. 0 ) cycle

        do j=1,sets(i)%ngeos
            ! ------------------------------------------------------------------
            if( .not. sets(i)%geo(j)%trg_probe_ene_loaded ) cycle

            nene = nene + 1

            err = 0.0d0
            select case(sets(i)%geo(j)%trg_probe_ene_mode)
                case(GEO_PROBE_ENE_REP)
                    emm = sets(i)%geo(j)%sapt_rep
                    err = emm - sets(i)%geo(j)%trg_probe_ene
                case(GEO_PROBE_ENE_TOT)
                    emm = sets(i)%geo(j)%sapt_ele + sets(i)%geo(j)%sapt_rep + sets(i)%geo(j)%sapt_dis
                    err = emm - sets(i)%geo(j)%trg_probe_ene
                case default
                    call ffdev_utils_exit(DEV_ERR,1,'Unsupported probe_mode in ffdev_err_probe_error!')
            end select
            serene = serene + sets(i)%geo(j)%weight * err**2

        end do
    end do

    if( nene .gt. 0 ) then
        error%probe_ene = sqrt(serene/real(nene))
    end if

end subroutine ffdev_err_probe_error

! ==============================================================================
! subroutine ffdev_err_probe_summary
! ==============================================================================

subroutine ffdev_err_probe_summary

    use ffdev_targetset_dat
    use ffdev_geometry
    use ffdev_err_probe_dat
    use ffdev_utils

    implicit none
    integer             :: i,j,nene
    logical             :: printsum
    real(DEVDP)         :: serene,errtot,emm,err
    ! --------------------------------------------------------------------------

    printsum = .false.
    do i=1,nsets
        do j=1,sets(i)%ngeos
            if( sets(i)%geo(j)%trg_probe_ene_loaded .and. ( sets(i)%top%probe_size .ne. 0) ) then
                printsum = .true.
            end if
        end do
    end do
    if( .not. printsum ) return

    write(DEV_OUT,*)
    write(DEV_OUT,5)
    write(DEV_OUT,10)
    write(DEV_OUT,20)

    serene = 0.0
    nene = 0

    do i=1,nsets
        printsum = .false.
        if( sets(i)%top%probe_size .eq. 0 ) cycle

        do j=1,sets(i)%ngeos
            if( .not. sets(i)%geo(j)%trg_probe_ene_loaded ) cycle

            printsum = .true.
            nene = nene + 1

            err = 0.0d0
            select case(sets(i)%geo(j)%trg_probe_ene_mode)
                case(GEO_PROBE_ENE_REP)
                    emm = sets(i)%geo(j)%sapt_rep
                    err = emm - sets(i)%geo(j)%trg_probe_ene
                case(GEO_PROBE_ENE_TOT)
                    emm = sets(i)%geo(j)%sapt_ele + sets(i)%geo(j)%sapt_rep + sets(i)%geo(j)%sapt_dis
                    err = emm - sets(i)%geo(j)%trg_probe_ene
                case default
                    call ffdev_utils_exit(DEV_ERR,1,'Unsupported probe_mode in ffdev_err_probe_summary!')
            end select
            serene = serene + sets(i)%geo(j)%weight * err**2

            write(DEV_OUT,30) i, j, sets(i)%geo(j)%weight, &
                              sets(i)%geo(j)%trg_probe_ene, emm, err
        end do
        if( printsum ) write(DEV_OUT,20)
    end do

    errtot = 0.0
    if( nene .gt. 0 ) then
        errtot  = sqrt(serene/real(nene))
    end if

    write(DEV_OUT,40)  errtot
    write(DEV_OUT,45)  ProbeErrorWeight*errtot

 5 format('# Probe errors')
10 format('# SET GeoID Weight   ELE(TGR)    ELE(MM) abs E(Err)')
20 format('# --- ----- ------ ---------- ---------- ----------')
30 format(I5,1X,I5,1X,F6.3,1X,F10.3,1X,F10.3,1X,F10.3)

40 format('# Final error (weighted per geometry)  = ',1X,F10.3)
45 format('# Final error (all weights)            = ',1X,F10.3)

end subroutine ffdev_err_probe_summary

! ------------------------------------------------------------------------------

end module ffdev_err_probe

