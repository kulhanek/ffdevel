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

module ffdev_errors

use ffdev_errors_dat
use ffdev_constants

contains

! ==============================================================================
! subroutine ffdev_errors_init_all
! ==============================================================================

subroutine ffdev_errors_init_all()

    use ffdev_err_bonds
    use ffdev_err_angles
    use ffdev_err_dihedrals
    use ffdev_err_impropers
    use ffdev_err_nbdists
    use ffdev_err_energy
    use ffdev_err_rmsd
    use ffdev_err_ihess
    use ffdev_errors_dat
    use ffdev_utils

    implicit none
    ! --------------------------------------------------------------------------

    ! clear what should be calculated
    errors_calc_hess = .false.
    errors_calc_grad = .false.
    errors_calc_ene  = .false.

    ! reset all setup
    call ffdev_err_energy_init()
    call ffdev_err_bonds_init()
    call ffdev_err_angles_init()
    call ffdev_err_dihedrals_init()
    call ffdev_err_impropers_init()
    call ffdev_err_nbdists_init()
    call ffdev_err_rmsd_init()

end subroutine ffdev_errors_init_all

! ==============================================================================
! subroutine ffdev_errors_error_only
! ==============================================================================

subroutine ffdev_errors_error_only(error)

    use ffdev_err_bonds_dat
    use ffdev_err_bonds

    use ffdev_err_angles_dat
    use ffdev_err_angles

    use ffdev_err_dihedrals_dat
    use ffdev_err_dihedrals

    use ffdev_err_nbdists_dat
    use ffdev_err_nbdists

    use ffdev_err_energy_dat
    use ffdev_err_energy

    use ffdev_err_ihess_dat
    use ffdev_err_ihess

    use ffdev_err_impropers_dat
    use ffdev_err_impropers

    use ffdev_err_rmsd_dat
    use ffdev_err_rmsd

    use ffdev_timers

    implicit none
    type(FFERROR_TYPE)  :: error
    ! --------------------------------------------------------------------------

    call ffdev_timers_start_timer(FFDEV_ERRORS_TIMER)

    error%total = 0.0d0
    error%energy = 0.0d0
    error%bonds = 0.0d0
    error%angles = 0.0d0
    error%dihedrals = 0.0d0
    error%impropers = 0.0d0
    error%nbdists = 0.0d0
    error%rmsd = 0.0d0

! energy based errors
    if( EnableEnergyError ) then
        call ffdev_err_energy_error(error)
        error%total = error%total + error%energy*EnergyErrorWeight
    end if

! geometry based errors
    if( EnableBondsError ) then
        call ffdev_err_bonds_error(error)
        error%total = error%total + error%bonds*BondErrorsWeight
    end if

    if( EnableAnglesError ) then
        call ffdev_err_angles_error(error)
        error%total = error%total + error%angles*AngleErrorsWeight
    end if

    if( EnableDihedralsError ) then
        call ffdev_err_dihedrals_error(error)
        error%total = error%total + error%dihedrals*DihedralsErrorWeight
    end if

    if( EnableImpropersError ) then
        call ffdev_err_impropers_error(error)
        error%total = error%total + error%impropers*ImpropersErrorWeight
    end if

    if( EnableNBDistsError ) then
        call ffdev_err_nbdists_error(error)
        error%total = error%total + error%nbdists*NBDistsErrorWeight
    end if

    if( EnableRMSDError ) then
        call ffdev_err_rmsd_error(error)
        error%total = error%total + error%rmsd*RMSDErrorWeight
    end if

    call ffdev_timers_stop_timer(FFDEV_ERRORS_TIMER)

end subroutine ffdev_errors_error_only

!===============================================================================
! subroutine ffdev_errors_ffopt_header_I
!===============================================================================

subroutine ffdev_errors_ffopt_header_I()

    use ffdev_err_bonds_dat
    use ffdev_err_angles_dat
    use ffdev_err_dihedrals_dat
    use ffdev_err_nbdists_dat
    use ffdev_err_energy_dat
    use ffdev_err_ihess_dat
    use ffdev_err_impropers_dat
    use ffdev_err_rmsd_dat

    implicit none
    ! --------------------------------------------------------------------------

    if( EnableEnergyError ) then
        write(DEV_OUT,30,ADVANCE='NO')
    end if
    if( EnableBondsError ) then
        write(DEV_OUT,33,ADVANCE='NO')
    end if
    if( EnableAnglesError ) then
        write(DEV_OUT,34,ADVANCE='NO')
    end if
    if( EnableDihedralsError ) then
        write(DEV_OUT,35,ADVANCE='NO')
    end if
    if( EnableImpropersError ) then
        write(DEV_OUT,38,ADVANCE='NO')
    end if
    if( EnableNBDistsError ) then
        write(DEV_OUT,36,ADVANCE='NO')
    end if
    if( EnableRMSDError ) then
        write(DEV_OUT,39,ADVANCE='NO')
    end if

 30 format('       Energy')
 33 format('        Bonds')
 34 format('       Angles')
 35 format('    Dihedrals')
 36 format('       d(NBs)')
 38 format('    Impropers')
 39 format('         RMSD')

end subroutine ffdev_errors_ffopt_header_I

!===============================================================================
! subroutine ffdev_errors_ffopt_header_II
!===============================================================================

subroutine ffdev_errors_ffopt_header_II()

    use ffdev_err_bonds_dat
    use ffdev_err_angles_dat
    use ffdev_err_dihedrals_dat
    use ffdev_err_nbdists_dat
    use ffdev_err_energy_dat
    use ffdev_err_impropers_dat
    use ffdev_err_rmsd_dat

    implicit none
    ! --------------------------------------------------------------------------

    if( EnableEnergyError ) then
        write(DEV_OUT,50,ADVANCE='NO')
    end if
    if( EnableBondsError ) then
        write(DEV_OUT,50,ADVANCE='NO')
    end if
    if( EnableAnglesError ) then
        write(DEV_OUT,50,ADVANCE='NO')
    end if
    if( EnableDihedralsError ) then
        write(DEV_OUT,50,ADVANCE='NO')
    end if
    if( EnableImpropersError ) then
        write(DEV_OUT,50,ADVANCE='NO')
    end if
    if( EnableNBDistsError ) then
        write(DEV_OUT,50,ADVANCE='NO')
    end if
    if( EnableRMSDError ) then
        write(DEV_OUT,50,ADVANCE='NO')
    end if

 50 format(' ------------')

end subroutine ffdev_errors_ffopt_header_II

!===============================================================================
! subroutine ffdev_errors_ffopt_header_II
!===============================================================================

subroutine ffdev_errors_ffopt_results(error)

    use ffdev_err_bonds_dat
    use ffdev_err_angles_dat
    use ffdev_err_dihedrals_dat
    use ffdev_err_nbdists_dat
    use ffdev_err_energy_dat
    use ffdev_err_impropers_dat
    use ffdev_err_rmsd_dat

    implicit none
    type(FFERROR_TYPE)  :: error
    ! -----------------------------------------------------------------------------

    if( EnableEnergyError ) then
        write(DEV_OUT,15,ADVANCE='NO') error%energy
    end if
    if( EnableBondsError ) then
        write(DEV_OUT,15,ADVANCE='NO') error%bonds
    end if
    if( EnableAnglesError ) then
        write(DEV_OUT,15,ADVANCE='NO') error%angles
    end if
    if( EnableDihedralsError ) then
        write(DEV_OUT,15,ADVANCE='NO') error%dihedrals
    end if
    if( EnableImpropersError ) then
        write(DEV_OUT,15,ADVANCE='NO') error%impropers
    end if
    if( EnableNBDistsError ) then
        write(DEV_OUT,15,ADVANCE='NO') error%nbdists
    end if
    if( EnableRMSDError ) then
        write(DEV_OUT,15,ADVANCE='NO') error%rmsd
    end if

 15 format(1X,E12.5)

end subroutine ffdev_errors_ffopt_results

! ==============================================================================
! subroutine ffdev_errors_summary
! ==============================================================================

subroutine ffdev_errors_summary(final)

    use ffdev_targetset_dat
    use ffdev_geometry
    use ffdev_utils

    use ffdev_err_bonds_dat
    use ffdev_err_angles_dat
    use ffdev_err_dihedrals_dat
    use ffdev_err_nbdists_dat
    use ffdev_err_energy_dat
    use ffdev_err_impropers_dat
    use ffdev_err_rmsd_dat

    use ffdev_err_bonds
    use ffdev_err_angles
    use ffdev_err_dihedrals
    use ffdev_err_nbdists
    use ffdev_err_energy
    use ffdev_err_impropers
    use ffdev_err_rmsd

    implicit none
    logical     :: final, anyrep, printme, printsum
    integer     :: i,j
    ! --------------------------------------------------------------------------

    if( .not. (PrintEnergyErrorSummary .or. &
            PrintBondsErrorSummary .or. PrintAnglesErrorSummary .or. PrintDihedralsErrorSummary .or. &
            PrintImpropersErrorSummary .or. &
            PrintNBDistsErrorSummary .or. PrintRMSDErrorSummary) ) then
        ! no error to report
        return
    end if

    write(DEV_OUT,*)
    write(DEV_OUT,1)
    if( final ) then
        call ffdev_utils_heading(DEV_OUT,'Final Error Statistics',':')
    else
        call ffdev_utils_heading(DEV_OUT,'Initial Error Statistics',':')
    end if
    write(DEV_OUT,1)

    ! summary pers sets
    if( PrintEnergyErrorSummary ) then
        write(DEV_OUT,*)
        write(DEV_OUT,10)

        do i=1,nsets
            printme = .false.
            if( PrintEnergyErrorSummary ) then
                printsum = .false.
                call ffdev_err_energy_summary(sets(i),printsum)
                printme = printme .or. printsum
            end if
            if( PrintRMSDErrorSummary ) then
                printsum = .false.
                call ffdev_err_rmsd_summary(sets(i),printsum)
                printme = printme .or. printsum
            end if

            if( .not. printme ) cycle

            printsum = .true.

            write(DEV_OUT,*)
            write(DEV_OUT,5) i
            if( PrintEnergyErrorSummary ) then
                call ffdev_err_energy_summary(sets(i),printsum)
            end if
            if( PrintRMSDErrorSummary ) then
                call ffdev_err_rmsd_summary(sets(i),printsum)
            end if
        end do
    end if

    ! summary per points
    if( PrintBondsErrorSummary .or. PrintAnglesErrorSummary .or. PrintDihedralsErrorSummary .or. &
        PrintImpropersErrorSummary .or. &
        PrintNBDistsErrorSummary .or. PrintRMSDErrorSummary) then

        write(DEV_OUT,*)
        write(DEV_OUT,20)

        do i=1,nsets
            do j=1,sets(i)%ngeos
                printme = .false.
                if( PrintBondsErrorSummary ) then
                    printsum = .false.
                    call ffdev_err_bonds_summary(sets(i)%top,sets(i)%geo(j),printsum)
                    printme = printme .or. printsum
                end if
                if( PrintAnglesErrorSummary ) then
                    printsum = .false.
                    call ffdev_err_angles_summary(sets(i)%top,sets(i)%geo(j),printsum)
                    printme = printme .or. printsum
                end if
                if( PrintDihedralsErrorSummary ) then
                    printsum = .false.
                    call ffdev_err_dihedrals_summary(sets(i)%top,sets(i)%geo(j),printsum)
                    printme = printme .or. printsum
                end if
                if( PrintImpropersErrorSummary ) then
                    printsum = .false.
                    call ffdev_err_impropers_summary(sets(i)%top,sets(i)%geo(j),printsum)
                    printme = printme .or. printsum
                end if
                if( PrintNBDistsErrorSummary ) then
                    printsum = .false.
                    call ffdev_err_nbdists_summary(sets(i)%top,sets(i)%geo(j),printsum)
                    printme = printme .or. printsum
                end if

                if( .not. printme ) cycle

                write(DEV_OUT,*)
                write(DEV_OUT,6) i,j
                printsum = .true.
                if( PrintBondsErrorSummary ) then
                    call ffdev_err_bonds_summary(sets(i)%top,sets(i)%geo(j),printsum)
                end if
                if( PrintAnglesErrorSummary ) then
                    call ffdev_err_angles_summary(sets(i)%top,sets(i)%geo(j),printsum)
                end if
                if( PrintDihedralsErrorSummary ) then
                    call ffdev_err_dihedrals_summary(sets(i)%top,sets(i)%geo(j),printsum)
                end if
                if( PrintImpropersErrorSummary ) then
                    call ffdev_err_impropers_summary(sets(i)%top,sets(i)%geo(j),printsum)
                end if
                if( PrintNBDistsErrorSummary ) then
                    call ffdev_err_nbdists_summary(sets(i)%top,sets(i)%geo(j),printsum)
                end if
            end do
        end do
    end if

 1 format('# ==============================================================================')
 5 format('== [SET] #',I2.2,' ===================================================================')
 6 format('== [SET] #',I2.2,' / PTS',I6.6' =======================================================')
10 format('== # SUMMARY PER SETS #')
20 format('== # SUMMARY PER POINTS #')

end subroutine ffdev_errors_summary

! ------------------------------------------------------------------------------

end module ffdev_errors
