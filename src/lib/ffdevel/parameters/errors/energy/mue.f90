! ==============================================================================
! This file is part of FFDevel.
!    Copyright (C) 2026 Petr Kulhanek, kulhanek@chemi.muni.cz
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
use ffdev_errors_dat

!===============================================================================

! MUE - Maximum Unsigned Energy error

type, extends(ErrorFceType) :: TypeEFTMUE
    contains
        ! executive methods
        procedure   :: set_title_errfce         => ffdev_err_mue_set_title
        procedure   :: calc_errfce              => ffdev_err_mue_error
end type TypeEFTMUE

contains

!===============================================================================
! Subroutine:  ffdev_err_mue_set_title
!===============================================================================

subroutine ffdev_err_mue_set_title(err_item)

    implicit none
    class(TypeEFTMUE)   :: err_item
    ! --------------------------------------------------------------------------

    err_item%Title = 'MUE'

end subroutine ffdev_err_mue_set_title

! ==============================================================================
! subroutine ffdev_err_mue
! ==============================================================================

subroutine ffdev_err_mue_error(err_item,opterr)

    use ffdev_targetset
    use ffdev_targetset_dat
    use ffdev_utils
    use ffdev_geometry
    use ffdev_errors_dat

    implicit none
    class(TypeEFTMUE)   :: err_item
    logical             :: opterr
    ! --------------------------------------------
    integer             :: i,j
    ! --------------------------------------------------------------------------

    err_item%RepValue = 0.0d0
    err_item%OptValue = 0.0d0
    if( .not. err_item%Enabled ) then
        if( opterr ) return
    end if

    do i=1,nsets
        ! use only sets, which can provide reliable energy
        if( .not. ( (sets(i)%nrefs .ge. 1) .or. (sets(i)%top%probe_size .gt. 0) ) ) cycle

        do j=1,sets(i)%ngeos
            ! ------------------------------------------------------------------
            if( .not. sets(i)%geo(j)%trg_ene_loaded ) cycle

            if( abs(sets(i)%geo(j)%total_ene - sets(i)%geo(j)%trg_energy) .gt. err_item%RepValue ) then
                err_item%RepValue = abs(sets(i)%geo(j)%total_ene - sets(i)%geo(j)%trg_energy)
            end if

        end do
    end do

    if( err_item%ScaleFac .gt. 0 ) then
        err_item%OptValue = err_item%Weight * err_item%RepValue**2 / err_item%ScaleFac**2
    end if

end subroutine ffdev_err_mue_error

! ------------------------------------------------------------------------------

end module ffdev_err_mue


