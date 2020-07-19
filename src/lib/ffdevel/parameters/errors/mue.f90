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

module ffdev_err_mue

use ffdev_constants
use ffdev_variables

contains

! ==============================================================================
! subroutine ffdev_err_mue_init
! ==============================================================================

subroutine ffdev_err_mue_init

    use ffdev_err_mue_dat
    use ffdev_errors_dat

    implicit none
    ! --------------------------------------------------------------------------

    EnableMUEError       = .false.
    MUEErrorWeight       = 1.0

end subroutine ffdev_err_mue_init

! ==============================================================================
! subroutine ffdev_err_mue
! ==============================================================================

subroutine ffdev_err_mue_error(error)

    use ffdev_targetset
    use ffdev_targetset_dat
    use ffdev_utils
    use ffdev_geometry
    use ffdev_errors_dat
    use ffdev_err_mue_dat

    implicit none
    type(FFERROR_TYPE)  :: error
    ! --------------------------------------------
    integer             :: i,j
    real(DEVDP)         :: mue
    ! --------------------------------------------------------------------------

    error%mue = 0.0d0

    do i=1,nsets
        ! use only sets, which can provide reliable energy
        if( .not. ( (sets(i)%nrefs .ge. 1) .or. (sets(i)%top%probe_size .gt. 0) ) ) cycle

        do j=1,sets(i)%ngeos
            ! ------------------------------------------------------------------
            if( .not. sets(i)%geo(j)%trg_ene_loaded ) cycle

            if( abs(sets(i)%geo(j)%total_ene - sets(i)%geo(j)%trg_energy) .gt. error%mue ) then
                error%mue = abs(sets(i)%geo(j)%total_ene - sets(i)%geo(j)%trg_energy)
            end if

        end do
    end do

end subroutine ffdev_err_mue_error

! ------------------------------------------------------------------------------

end module ffdev_err_mue


