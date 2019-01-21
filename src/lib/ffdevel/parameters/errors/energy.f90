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

module ffdev_err_energy

use ffdev_constants

contains

! ==============================================================================
! subroutine ffdev_err_energy_error
! ==============================================================================

subroutine ffdev_err_energy_error(error)

    use ffdev_targetset
    use ffdev_targetset_dat
    use ffdev_utils
    use ffdev_geometry
    use ffdev_errors_dat

    implicit none
    type(FFERROR_TYPE)  :: error
    ! --------------------------------------------
    integer             :: i,j,nene
    real(DEVDP)         :: err,seterrene
    ! --------------------------------------------------------------------------

    error%energy = 0.0d0

    seterrene = 0.0
    nene = 0

    do i=1,nsets
        do j=1,sets(i)%ngeos
            ! ------------------------------------------------------------------
            if( sets(i)%geo(j)%trg_ene_loaded ) then
                nene = nene + 1
                err = sets(i)%geo(j)%total_ene - sets(i)%offset - sets(i)%geo(j)%trg_energy
                seterrene = seterrene + sets(i)%geo(j)%weight * err**2
            end if
        end do
    end do

    if( nene .gt. 0 ) then
        error%energy = sqrt(seterrene/real(nene))
    end if

end subroutine ffdev_err_energy_error

! ------------------------------------------------------------------------------

end module ffdev_err_energy


