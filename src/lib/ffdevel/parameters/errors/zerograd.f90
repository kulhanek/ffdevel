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

module ffdev_err_zerograd

use ffdev_constants
use ffdev_variables

contains

! ==============================================================================
! subroutine ffdev_err_zerograd_init
! ==============================================================================

subroutine ffdev_err_zerograd_init

    use ffdev_err_zerograd_dat
    use ffdev_errors_dat

    implicit none
    ! --------------------------------------------------------------------------

    EnableZeroGradError       = .false.
    PrintZeroGradErrorSummary = .false.
    ZeroGradErrorWeight       = 1.0

end subroutine ffdev_err_zerograd_init

! ==============================================================================
! subroutine ffdev_err_zerograd_error
! ==============================================================================

subroutine ffdev_err_zerograd_error(error)

    use ffdev_targetset
    use ffdev_targetset_dat
    use ffdev_utils
    use ffdev_geometry
    use ffdev_errors_dat
    use ffdev_err_zerograd_dat

    implicit none
    type(FFERROR_TYPE)  :: error
    ! --------------------------------------------
    integer             :: i,j,k
    real(DEVDP)         :: grms,totgrms,nele
    ! --------------------------------------------------------------------------

    error%zerograd = 0.0d0

    nele = 0
    totgrms = 0.0d0
    do i=1,nsets
        if( sets(i)%top%probe_size .ne. 0 ) cycle   ! skip probes
        if( sets(i)%top%natoms .le. 1 ) cycle       ! at least two atoms for gradient

        do j=1,sets(i)%ngeos
            grms = 0.0d0
            do k=1,sets(i)%top%natoms
                grms = grms + sets(i)%geo(j)%grd(1,k)**2 + sets(i)%geo(j)%grd(2,k)**2  + sets(i)%geo(j)%grd(3,k)**2
            end do
            if( sets(i)%top%natoms .gt. 0 ) then
                grms = sqrt(grms/real(3 * sets(i)%top%natoms))
            end if
            totgrms = totgrms + sets(i)%geo(j)%weight * grms ** 2
            nele = nele + 1
        end do
    end do

    if( nele.gt. 0 ) then
        error%zerograd = sqrt(totgrms/nele)
    end if

end subroutine ffdev_err_zerograd_error

! ==============================================================================
! subroutine ffdev_err_zerograd_summary
! ==============================================================================

subroutine ffdev_err_zerograd_summary

    use ffdev_targetset_dat
    use ffdev_geometry
    use ffdev_err_zerograd_dat

    implicit none
    real(DEVDP)         :: grms,totgrms,nele
    integer             :: i,j,k
    logical             :: printsum
    ! --------------------------------------------------------------------------

    printsum = .false.
    do i=1,nsets
        do j=1,sets(i)%ngeos
            printsum = .true.
        end do
    end do
    if( .not. printsum ) return

    write(DEV_OUT,*)
    write(DEV_OUT,5)
    write(DEV_OUT,10)
    write(DEV_OUT,20)

    totgrms = 0.0d0

    nele = 0
    totgrms = 0.0d0
    do i=1,nsets
        if( sets(i)%top%probe_size .ne. 0 ) cycle   ! skip probes
        if( sets(i)%top%natoms .le. 1 ) cycle       ! at least two atoms for gradient

        do j=1,sets(i)%ngeos
            printsum = .true.
            grms = 0.0d0
            do k=1,sets(i)%top%natoms
                grms = grms + sets(i)%geo(j)%grd(1,k)**2 + sets(i)%geo(j)%grd(2,k)**2  + sets(i)%geo(j)%grd(3,k)**2
            end do
            if( sets(i)%top%natoms .gt. 0 ) then
                grms = sqrt(grms/real(3 * sets(i)%top%natoms))
            end if
            totgrms = totgrms + sets(i)%geo(j)%weight * grms ** 2
            nele = nele + 1
            write(DEV_OUT,30) i, j, sets(i)%geo(j)%weight, grms
        end do
        if( printsum ) write(DEV_OUT,20)
    end do

    if( nele.gt. 0 ) then
        totgrms = sqrt(totgrms/nele)
    end if

    write(DEV_OUT,40)  totgrms
    write(DEV_OUT,45)  ZeroGradErrorWeight*totgrms

 5 format('# Zero gradient errors')
10 format('# SET GeoID Weight   GRMS(MM)')
20 format('# --- ----- ------ ----------')
30 format(I5,1X,I5,1X,F6.3,1X,F10.3)
40 format('# Final error (weighted per geometry) =  ',F10.3)
45 format('# Final error (all weights)           =  ',F10.3)

end subroutine ffdev_err_zerograd_summary

! ------------------------------------------------------------------------------

end module ffdev_err_zerograd


