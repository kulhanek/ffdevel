! ==============================================================================
! This file is part of FFDevel.
!    Copyright (C) 2026 Petr Kulhanek, kulhanek@chemi.muni.cz
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

module ffdev_errors

use ffdev_errors_dat
use ffdev_constants
use ffdev_variables

contains

! ==============================================================================
! subroutine ffdev_errors_init
! ==============================================================================

subroutine ffdev_errors_init()

    use ffdev_errors_dat

    implicit none
    integer             :: i
    ! --------------------------------------------------------------------------

    ! clear what should be calculated
    errors_calc_ene     = .false.
    errors_calc_sapt    = .false.
    errors_calc_grad    = .false.
    errors_calc_hess    = .false.

    ! init all error
    do i=1,NumOfErrorFces
        call ErrorFceList(i)%ErrFce%init_errfce()
    end do

end subroutine ffdev_errors_init

! ==============================================================================
! subroutine ffdev_errors_error_setup_domains
! opterr - calculate minimum for error optimization
! ==============================================================================

subroutine ffdev_errors_error_setup_domains(opterr)

    use ffdev_errors_dat

    implicit none
    logical         :: opterr
    ! --------------------------------------------
    integer         :: i
    ! --------------------------------------------------------------------------

    ! clear what should be calculated
    errors_calc_ene     = .false.
    errors_calc_sapt    = .false.
    errors_calc_grad    = .false.
    errors_calc_hess    = .false.

    ! setup domains
    do i=1,NumOfErrorFces
        call ErrorFceList(i)%ErrFce%setup_domains_errfce(opterr)
    end do

! to be moved to individual classes
!    errors_calc_ene  = EnableEnergyError .or. EnableProbeError
!    errors_calc_sapt = EnableSAPTError
!    errors_calc_grad = EnableZeroGradError
!
!    if( .not. opterr ) then
!        errors_calc_ene  = errors_calc_ene  .or. PrintEnergyErrorSummary .or. PrintProbeErrorSummary
!        errors_calc_sapt = errors_calc_sapt .or. PrintSAPTErrorSummary
!        errors_calc_grad = errors_calc_grad .or. PrintZeroGradErrorSummary
!    end if

end subroutine ffdev_errors_error_setup_domains

! ==============================================================================
! subroutine ffdev_errors_error_only
! opterr - calculate minimum for error optimization
! ==============================================================================

subroutine ffdev_errors_error_only(errfcetot,opterr)

    use ffdev_errors_dat
    use ffdev_timers

    implicit none
    real(DEVDP)     :: errfcetot
    logical         :: opterr
    ! --------------------------------------------
    integer         :: i
    ! --------------------------------------------------------------------------

    call ffdev_timers_start_timer(FFDEV_ERRORS_TIMER)

    errfcetot = 0.0d0

    ! get individual errors
    do i=1,NumOfErrorFces
        call ErrorFceList(i)%ErrFce%calc_errfce(opterr)
        errfcetot = errfcetot + ErrorFceList(i)%ErrFce%OptValue
    end do

    call ffdev_timers_stop_timer(FFDEV_ERRORS_TIMER)

end subroutine ffdev_errors_error_only

!===============================================================================
! subroutine ffdev_errors_ffopt_header_I
!===============================================================================

subroutine ffdev_errors_ffopt_header_I()

    use ffdev_errors_dat

    implicit none
    integer             :: i
    ! --------------------------------------------------------------------------

    ! get individual errors
    do i=1,NumOfErrorFces
        write(DEV_OUT,10,ADVANCE='NO') trim(ErrorFceList(i)%ErrFce%Title)
    end do

 10 format(1X,A12)

end subroutine ffdev_errors_ffopt_header_I

!===============================================================================
! subroutine ffdev_errors_ffopt_header_II
!===============================================================================

subroutine ffdev_errors_ffopt_header_II()

    use ffdev_errors_dat

    implicit none
    integer             :: i
    ! --------------------------------------------------------------------------

    do i=1,NumOfErrorFces
        write(DEV_OUT,10,ADVANCE='NO')
    end do

 10 format(' ------------')

end subroutine ffdev_errors_ffopt_header_II

!===============================================================================
! subroutine ffdev_errors_ffopt_header_scale_fac
!===============================================================================

subroutine ffdev_errors_ffopt_header_scale_fac()

    use ffdev_errors_dat

    implicit none
    integer             :: i
    ! --------------------------------------------------------------------------

    do i=1,NumOfErrorFces
        write(DEV_OUT,10,ADVANCE='NO') ErrorFceList(i)%ErrFce%ScaleFac
    end do

10 format(1X,E12.5)

end subroutine ffdev_errors_ffopt_header_scale_fac

!===============================================================================
! subroutine ffdev_errors_ffopt_header_weight
!===============================================================================

subroutine ffdev_errors_ffopt_header_weight()

    use ffdev_errors_dat

    implicit none
    integer             :: i
    ! --------------------------------------------------------------------------

    do i=1,NumOfErrorFces
        write(DEV_OUT,10,ADVANCE='NO') ErrorFceList(i)%ErrFce%Weight
    end do

10 format(1X,E12.5)

end subroutine ffdev_errors_ffopt_header_weight

!===============================================================================
! subroutine ffdev_errors_ffopt_results
!===============================================================================

subroutine ffdev_errors_ffopt_results()

    use ffdev_errors_dat

    implicit none
    integer             :: i
    ! --------------------------------------------------------------------------

    do i=1,NumOfErrorFces
        write(DEV_OUT,10,ADVANCE='NO') ErrorFceList(i)%ErrFce%RepValue
    end do

10 format(1X,E12.5)

end subroutine ffdev_errors_ffopt_results

! ==============================================================================
! subroutine ffdev_errors_summary
! ==============================================================================

subroutine ffdev_errors_summary(logmode)

    use ffdev_errors_dat
    use ffdev_targetset_dat
    use ffdev_utils

    implicit none
    integer     :: logmode
    ! --------------------------------------------
    logical     :: printflag, printsum
    integer     :: i,j,k
    ! --------------------------------------------------------------------------

    printflag = .false.
    do k=1,NumOfErrorFces
        printflag = printflag .or. ErrorFceList(k)%ErrFce%PrintSummary
    end do

    if( .not. printflag ) then
        ! no error to report
        return
    end if

    write(DEV_OUT,*)
    write(DEV_OUT,1)
    select case(logmode)
        case(SMMLOG_INITIAL)
            call ffdev_utils_heading(DEV_OUT,'Initial Error Statistics',':')
        case(SMMLOG_INTERMEDIATE)
            call ffdev_utils_heading(DEV_OUT,'Intermediate Error Statistics',':')
        case(SMMLOG_FINAL)
            call ffdev_utils_heading(DEV_OUT,'Final Error Statistics',':')
    end select
    write(DEV_OUT,1)

    ! individual summaries
    do k=1,NumOfErrorFces
        call ErrorFceList(k)%ErrFce%print_individual_summary_errfce()
    end do

    ! summary per sets
    printflag = .false.
    do i=1,nsets
        do k=1,NumOfErrorFces
            printsum = .false.
            call ErrorFceList(k)%ErrFce%print_set_summary_errfce(sets(i),printsum)
            printflag = printflag .or. printsum
        end do
    end do

    if( printflag ) then
        write(DEV_OUT,*)
        write(DEV_OUT,10)

        do i=1,nsets
            do k=1,NumOfErrorFces
                printsum = .false.
                call ErrorFceList(k)%ErrFce%print_set_summary_errfce(sets(i),printsum)
                printflag = printflag .or. printsum

                if( .not. printflag ) cycle

                write(DEV_OUT,*)
                write(DEV_OUT,5) i
                printsum = .true.
                call ErrorFceList(k)%ErrFce%print_set_summary_errfce(sets(i),printsum)
            end do
        end do
    end if

    ! summary per points
    printflag = .false.
    do i=1,nsets
        do j=1,sets(i)%ngeos
            do k=1,NumOfErrorFces
                printsum = .false.
                call ErrorFceList(k)%ErrFce%print_pts_summary_errfce(sets(i)%top,sets(i)%geo(j),printsum)
                printflag = printflag .or. printsum
            end do
        end do
    end do ! <- do we need to print anything?

    if( printflag ) then
        write(DEV_OUT,*)
        write(DEV_OUT,20)

        do i=1,nsets
            do j=1,sets(i)%ngeos

                printflag = .false.
                do k=1,NumOfErrorFces
                    printsum = .false.
                    call ErrorFceList(k)%ErrFce%print_pts_summary_errfce(sets(i)%top,sets(i)%geo(j),printsum)
                    printflag = printflag .or. printsum
                end do  ! <- do we need to print something for the PTS?

                if( .not. printflag ) cycle

                write(DEV_OUT,*)
                write(DEV_OUT,6) i,j

                do k=1,NumOfErrorFces
                    printsum = .true.
                    call ErrorFceList(k)%ErrFce%print_pts_summary_errfce(sets(i)%top,sets(i)%geo(j),printsum)
                end do
            end do
        end do
    end if

 1 format('# ==============================================================================')
 5 format('== [SET] #',I2.2,' ===================================================================')
 6 format('== [SET#',I5.5,']/[GEO#',I6.6,'] ====================================================')
10 format('== # SUMMARY PER SETS #')
20 format('== # SUMMARY PER POINTS #')

end subroutine ffdev_errors_summary

! ------------------------------------------------------------------------------

end module ffdev_errors
