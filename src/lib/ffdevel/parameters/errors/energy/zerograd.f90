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

module ffdev_err_zerograd

use ffdev_constants
use ffdev_variables
use ffdev_errors_dat

!===============================================================================

type, extends(ErrorFceType) :: TypeEFTZeroGrad
    contains
        ! executive methods
        procedure   :: set_title_errfce         => ffdev_err_zerograd_set_title
        procedure   :: calc_errfce              => ffdev_err_zerograd_error
        procedure   :: print_individual_summary_errfce => ffdev_err_zerograd_summary
end type TypeEFTZeroGrad

contains

!===============================================================================
! Subroutine:  ffdev_err_zerograd_set_title
!===============================================================================

subroutine ffdev_err_zerograd_set_title(err_item)

    implicit none
    class(TypeEFTZeroGrad)    :: err_item
    ! --------------------------------------------------------------------------

    err_item%Title = 'ZeroGrd'

end subroutine ffdev_err_zerograd_set_title

! ==============================================================================
! subroutine ffdev_err_zerograd_error
! ==============================================================================

subroutine ffdev_err_zerograd_error(err_item,opterr)

    use ffdev_targetset
    use ffdev_targetset_dat
    use ffdev_utils
    use ffdev_geometry

    implicit none
    class(TypeEFTZeroGrad)  :: err_item
    logical                 :: opterr
    ! --------------------------------------------
    integer             :: i,j,k
    real(DEVDP)         :: grms,totgrms,nele
    ! --------------------------------------------------------------------------

    err_item%RepValue = 0.0d0
    err_item%OptValue = 0.0d0
    if( .not. err_item%Enabled ) then
        if( opterr ) return
    end if

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
        err_item%RepValue = sqrt(totgrms/nele)
    end if

    if( err_item%ScaleFac .gt. 0 ) then
        err_item%OptValue = err_item%Weight * err_item%RepValue**2 / err_item%ScaleFac**2
    end if

end subroutine ffdev_err_zerograd_error

! ==============================================================================
! subroutine ffdev_err_zerograd_summary
! ==============================================================================

subroutine ffdev_err_zerograd_summary(err_item)

    use ffdev_targetset_dat
    use ffdev_geometry

    implicit none
    class(TypeEFTZeroGrad)  :: err_item
    ! --------------------------------------------
    real(DEVDP)         :: grms,totgrms,totw
    integer             :: i,j,k
    logical             :: printsum
    ! --------------------------------------------------------------------------

    if( .not. err_item%PrintSummary ) return

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

    totw = 0
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
            totw = totw + sets(i)%geo(j)%weight
            write(DEV_OUT,30) i, j, sets(i)%geo(j)%weight, grms
        end do
        if( printsum ) write(DEV_OUT,20)
    end do

    if( totw .gt. 0 ) then
        totgrms = sqrt(totgrms/totw)
    end if

    write(DEV_OUT,40)  totgrms
    write(DEV_OUT,45)  err_item%Weight*totgrms

 5 format('# Zero gradient errors')
10 format('# SET  GeoID Weight   GRMS(MM)')
20 format('# --- ------ ------ ----------')
30 format(I5,1X,I6,1X,F6.3,1X,F10.3)
40 format('# Final error (weighted per geometry) =  ',F10.3)
45 format('# Final error (all weights)           =  ',F10.3)

end subroutine ffdev_err_zerograd_summary

! ------------------------------------------------------------------------------

end module ffdev_err_zerograd


