! ==============================================================================
! This file is part of FFDevel.
!    Copyright (C) 2019 Petr Kulhanek, kulhanek@chemi.muni.cz
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

module ffdev_err_nbdists

use ffdev_constants

contains

! ==============================================================================
! subroutine ffdev_err_nbdists_init
! ==============================================================================

subroutine ffdev_err_nbdists_init

    use ffdev_err_nbdists_dat

    implicit none
    ! --------------------------------------------------------------------------

    EnableNBDistanceError        = .false.
    PrintNBDistanceErrorSummary  = .false.
    NBDistanceErrorWeight        = 1.0

    NBDistanceSWPosition         = 4.0
    NBDistanceSWAlpha            = 1.0

end subroutine ffdev_err_nbdists_init

! ==============================================================================
! subroutine ffdev_err_nbdists_error
! ==============================================================================

subroutine ffdev_err_nbdists_error(error)

    use ffdev_targetset_dat
    use ffdev_utils   
    use ffdev_geometry
    use ffdev_errors_dat
    use ffdev_err_nbdists_dat

    implicit none
    type(FFERROR_TYPE)  :: error
    ! --------------------------------------------
    integer             :: i,j,q,ai,aj
    real(DEVDP)         :: err,seterrnbdists
    real(DEVDP)         :: d0,dt,sw,swsum
    ! --------------------------------------------------------------------------

    error%nbdists = 0.0

    seterrnbdists = 0.0
    swsum = 0

    do i=1,nsets
        do j=1,sets(i)%ngeos
            if( sets(i)%geo(j)%trg_crd_optimized .and. sets(i)%top%nfragments .gt. 1) then
                do q=1,sets(i)%top%nb_size
                    ai = sets(i)%top%nb_list(q)%ai
                    aj = sets(i)%top%nb_list(q)%aj

                    if( sets(i)%top%atoms(ai)%frgid .eq. sets(i)%top%atoms(aj)%frgid ) cycle

                    d0 = ffdev_geometry_get_length(sets(i)%geo(j)%crd,ai,aj)
                    dt = ffdev_geometry_get_length(sets(i)%geo(j)%trg_crd,ai,aj)
                    err = d0 - dt

                    ! calculate switch function
                    sw = 1.0d0 / (1.0d0 + exp( NBDistanceSWAlpha*(dt - NBDistanceSWPosition) ) )
                    swsum = swsum + sw
                    err = err * sw
                    seterrnbdists = seterrnbdists + sets(i)%geo(j)%weight * err**2
                end do
            end if
        end do
    end do

    if( swsum .gt. 0 ) then
        error%nbdists = sqrt(seterrnbdists/real(swsum))
    end if 

end subroutine ffdev_err_nbdists_error

! ------------------------------------------------------------------------------

end module ffdev_err_nbdists
