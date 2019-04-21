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

module ffdev_err_freqs

use ffdev_constants

contains

! ==============================================================================
! subroutine ffdev_err_freqs_init
! ==============================================================================

subroutine ffdev_err_freqs_init

    use ffdev_err_freqs_dat

    implicit none
    ! --------------------------------------------------------------------------

    EnableFreqsError        = .false.
    DebugFreqError          = .false.
    PrintFreqsErrorSummary  = .false.
    FreqsErrorWeight        = 1.0d0
    FreqsMaxAngle           = 40.0d0
    FreqsMaxRMSD            = 0.5d0
    FreqsErrorMode          = FREQS_SUPERIMPOSE_GEO

end subroutine ffdev_err_freqs_init

! ==============================================================================
! subroutine ffdev_err_freqs_error
! ==============================================================================

subroutine ffdev_err_freqs_error(error)

    use ffdev_targetset
    use ffdev_targetset_dat
    use ffdev_utils
    use ffdev_geometry
    use ffdev_errors_dat
    use ffdev_err_freqs_dat
    use ffdev_hessian_utils

    implicit none
    type(FFERROR_TYPE)  :: error
    ! --------------------------------------------
    integer             :: i,j,q,nfreqs
    real(DEVDP)         :: err,seterrfreqs,diff
    real(DEVDP)         :: rmsd
    ! --------------------------------------------------------------------------

    error%freqs = 0.0d0

    ! calculate error
    seterrfreqs = 0.0
    nfreqs = 0

    do i=1,nsets
        do j=1,sets(i)%ngeos
            ! ------------------------------------------------------------------
            if( sets(i)%geo(j)%trg_crd_optimized .and. sets(i)%geo(j)%trg_hess_loaded ) then

                if( DebugFreqError ) then
                    write(DEV_OUT,*)
                    write(DEV_OUT,10) i,j
                end if

                select case(FreqsErrorMode)
                    case(FREQS_SUPERIMPOSE_GEO)
                        ! superimpose and then find mapping
                        call ffdev_hessian_superimpose_by_geo(sets(i)%geo(j))

                    case(FREQS_SUPERIMPOSE_MODE)
                        ! superimpose and then find mapping
                        call ffdev_hessian_superimpose_by_modes(sets(i)%geo(j))

                    case default
                        call ffdev_utils_exit(DEV_OUT,1,'Unsupported mode in ffdev_err_freqs_error!')
                end select

                ! calculate error
                do q=7,3*sets(i)%geo(j)%natoms
                    ! rmsd check
                    if( sets(i)%geo(j)%freq_t2s_rmsd(q) .gt. FreqsMaxRMSD ) cycle

                    ! angle check
                    if( sets(i)%geo(j)%freq_t2s_angles(q)*DEV_R2D .gt. FreqsMaxAngle ) cycle

                    err = sets(i)%geo(j)%freq(sets(i)%geo(j)%freq_t2s_map(q))-sets(i)%geo(j)%trg_freq(q)
                    seterrfreqs = seterrfreqs + sets(i)%geo(j)%weight * err**2
                    nfreqs = nfreqs + 1
                end do

                if( DebugFreqError ) then
                    call ffdev_hessian_print_mapping(.false.,sets(i)%geo(j)%natoms,sets(i)%geo(j)%trg_freq, &
                                                     sets(i)%geo(j)%trg_nmodes, &
                                                     sets(i)%geo(j)%freq,sets(i)%geo(j)%nmodes, &
                                                     sets(i)%geo(j)%freq_t2s_map,sets(i)%geo(j)%freq_t2s_angles, &
                                                     sets(i)%geo(j)%freq_t2s_rmsd)
                end if


            end if
        end do
    end do

    if( DebugFreqError ) then
        write(DEV_OUT,*)
    end if

    if( nfreqs .gt. 0 ) then
        error%freqs = sqrt(seterrfreqs/real(nfreqs))
    end if

    10 format('=== #',I2.2,'/',I3.3)
    20 format('    RMSD = ',F10.3)

end subroutine ffdev_err_freqs_error

! ==============================================================================
! subroutine ffdev_err_freqs_summary
! ==============================================================================

subroutine ffdev_err_freqs_summary(geo,printsum)

    use ffdev_geometry_dat
    use ffdev_hessian_utils

    implicit none
    type(GEOMETRY)  :: geo
    logical         :: printsum
    ! --------------------------------------------

    if( .not. geo%trg_freq_loaded ) return
    if( .not. geo%trg_crd_optimized ) return

    if( printsum .eqv. .false. ) then
        printsum = .true.
        return
    end if

    call ffdev_hessian_print_mapping(.false.,geo%natoms,geo%trg_freq,geo%trg_nmodes,&
                                     geo%freq,geo%nmodes,geo%freq_t2s_map,geo%freq_t2s_angles,&
                                     geo%freq_t2s_rmsd)

end subroutine ffdev_err_freqs_summary

! ------------------------------------------------------------------------------

end module ffdev_err_freqs
