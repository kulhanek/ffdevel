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

    EnableFreqsError      = .false.
    DebugFreqError       = .false.
    PrintFreqsErrorSummary= .false.
    FreqsErrorWeight      = 1.0d0
    FreqMaxNmodeAngle    = 50.0d0
    FreqMaxRMSD          = 5.0d0

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
    integer             :: i,j,q,nfreqs,k,l,i1,j1,i2,j2
    real(DEVDP)         :: err,seterrfreqs,diff,nmangle
    real(DEVDP)         :: f0,ft,rmsd
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

                ! superimpose and then find mapping
                call ffdev_hessian_superimpose_freqs(sets(i)%geo(j)%natoms,sets(i)%geo(j)%z,sets(i)%geo(j)%trg_crd, &
                                                     sets(i)%geo(j)%crd,sets(i)%geo(j)%nmodes,.false.,rmsd)

                if( DebugFreqError ) then
                    write(DEV_OUT,20) rmsd
                end if

                if( rmsd .gt. FreqMaxRMSD ) cycle ! rmsd is not OK

                call ffdev_hessian_find_mapping(sets(i)%geo(j)%natoms,sets(i)%geo(j)%trg_nmodes, &
                                                sets(i)%geo(j)%nmodes,sets(i)%geo(j)%freq_t2s_map)

                ! write(*,*) rmsd

                do q=7,3*sets(i)%geo(j)%natoms
                    err = sets(i)%geo(j)%freq(sets(i)%geo(j)%freq_t2s_map(q))-sets(i)%geo(j)%trg_freq(q)

                    i1 = (q-1) / 3 + 1
                    j1 = mod(q-1,3) + 1

                    i2 = (sets(i)%geo(j)%freq_t2s_map(q)-1) / 3 + 1
                    j2 = mod(sets(i)%geo(j)%freq_t2s_map(q)-1,3) + 1

                    nmangle = 0.0
                    do k=1,sets(i)%geo(j)%natoms
                        do l=1,3
                            nmangle = nmangle + sets(i)%geo(j)%trg_nmodes(l,k,j1,i1)*sets(i)%geo(j)%nmodes(l,k,j2,i2)
                        end do
                    end do
                    nmangle = acos(nmangle)*DEV_R2D

                    ! write(*,*) nmangle,err,i1,j1,i2,j2

                    if( (nmangle .lt. FreqMaxNmodeAngle) .or. ((180d0 - nmangle) .lt. FreqMaxNmodeAngle) ) then
                        seterrfreqs = seterrfreqs + sets(i)%geo(j)%weight * err**2
                        nfreqs = nfreqs + 1
                    end if
                end do
                if( DebugFreqError ) then
                    call ffdev_hessian_print_mapping(.false.,sets(i)%geo(j)%natoms,sets(i)%geo(j)%trg_freq, &
                                                     sets(i)%geo(j)%trg_nmodes, &
                                                     sets(i)%geo(j)%freq,sets(i)%geo(j)%nmodes, &
                                                     sets(i)%geo(j)%freq_t2s_map)
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

subroutine ffdev_err_freqs_summary(geo)

    use ffdev_geometry_dat
    use ffdev_hessian_utils

    implicit none
    type(GEOMETRY)     :: geo
    ! --------------------------------------------

    call ffdev_hessian_print_mapping(.false.,geo%natoms,geo%trg_freq,geo%trg_nmodes,&
                                     geo%freq,geo%nmodes,geo%freq_t2s_map)

end subroutine ffdev_err_freqs_summary

! ------------------------------------------------------------------------------

end module ffdev_err_freqs
