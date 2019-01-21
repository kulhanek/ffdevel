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

module ffdev_err_bonds

use ffdev_constants

contains

! ==============================================================================
! subroutine ffdev_err_bonds_error
! ==============================================================================

subroutine ffdev_err_bonds_error(error)

    use ffdev_targetset
    use ffdev_targetset_dat
    use ffdev_utils
    use ffdev_geometry
    use ffdev_errors_dat

    implicit none
    type(FFERROR_TYPE)  :: error
    ! --------------------------------------------
    integer             :: i,j,q,nbonds,ai,aj
    real(DEVDP)         :: err,seterrbonds
    real(DEVDP)         :: d0,dt
    ! --------------------------------------------------------------------------

    error%bonds = 0.0d0

    ! calculate error
    seterrbonds = 0.0
    nbonds = 0

    do i=1,nsets
        do j=1,sets(i)%ngeos
            ! ------------------------------------------------------------------
            if( sets(i)%geo(j)%trg_crd_optimized ) then
                do q=1,sets(i)%top%nbonds
                    ai = sets(i)%top%bonds(q)%ai
                    aj = sets(i)%top%bonds(q)%aj
                    d0 = ffdev_geometry_get_length(sets(i)%geo(j)%crd,ai,aj)
                    dt = ffdev_geometry_get_length(sets(i)%geo(j)%trg_crd,ai,aj)
                    nbonds = nbonds + 1
                    err = d0 - dt
                    seterrbonds = seterrbonds + sets(i)%geo(j)%weight * err**2
                end do
            end if
        end do
    end do

    ! geometry
    if( nbonds .gt. 0 ) then
        error%bonds = sqrt(seterrbonds/real(nbonds))
    end if

end subroutine ffdev_err_bonds_error

! ------------------------------------------------------------------------------

end module ffdev_err_bonds
