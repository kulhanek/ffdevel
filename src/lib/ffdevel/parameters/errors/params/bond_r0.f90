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

module ffdev_err_bond_r0

use ffdev_constants
use ffdev_variables
use ffdev_errors_dat

!===============================================================================

type, extends(ErrorFceType) :: TypeEFTBondR0
    contains
        ! executive methods
        procedure   :: set_title_errfce         => ffdev_err_bond_r0_set_title
        procedure   :: calc_errfce              => ffdev_err_bond_r0_error
end type TypeEFTBondR0

contains

!===============================================================================
! Subroutine:  ffdev_err_bond_r0_set_title
!===============================================================================

subroutine ffdev_err_bond_r0_set_title(err_item)

    implicit none
    class(TypeEFTBondR0)    :: err_item
    ! --------------------------------------------------------------------------

    err_item%Title = 'bond_r0'

end subroutine ffdev_err_bond_r0_set_title

! ==============================================================================
! subroutine ffdev_err_bond_r0_error
! opterr - calculate minimum for error optimization
! ==============================================================================

subroutine ffdev_err_bond_r0_error(err_item,opterr)

    use ffdev_targetset
    use ffdev_targetset_dat
    use ffdev_utils
    use ffdev_geometry
    use ffdev_errors_dat

    implicit none
    class(TypeEFTBondR0) :: err_item
    logical             :: opterr
    ! --------------------------------------------
    integer             :: i,j,q,ai,aj
    real(DEVDP)         :: err,seterrbonds,totw
    real(DEVDP)         :: d0,dt
    ! --------------------------------------------------------------------------

    err_item%RepValue = 0.0d0
    err_item%OptValue = 0.0d0
    if( .not. err_item%Enabled ) then
        if( opterr ) return
    end if

    ! calculate error
    seterrbonds = 0.0
    totw = 0

    do i=1,nsets
        do q=1,sets(i)%top%nbonds

            if( err_item%OnlyFFOpt ) then
                if( .not. sets(i)%top%bond_types(sets(i)%top%bonds(q)%bt)%ffoptactive ) cycle
            end if
            ai = sets(i)%top%bonds(q)%ai
            aj = sets(i)%top%bonds(q)%aj

            d0 = sets(i)%top%bond_types(sets(i)%top%bonds(q)%bt)%d0

            do j=1,sets(i)%ngeos
                if( .not. sets(i)%geo(j)%trg_crd_loaded ) cycle

                dt = ffdev_geometry_get_length(sets(i)%geo(j)%trg_crd,ai,aj)

                ! write(*,*) d0, dt
                err = d0 - dt
                seterrbonds = seterrbonds + sets(i)%geo(j)%weight * err**2
                totw = totw + sets(i)%geo(j)%weight
            end do
        end do
    end do

    ! geometry
    if( totw .gt. 0 ) then
        err_item%RepValue = sqrt(seterrbonds/totw)
    end if

    if( err_item%ScaleFac .gt. 0 ) then
        err_item%OptValue = err_item%Weight * err_item%RepValue**2 / err_item%ScaleFac**2
    end if

end subroutine ffdev_err_bond_r0_error

! ------------------------------------------------------------------------------

end module ffdev_err_bond_r0
