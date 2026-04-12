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

module ffdev_err_bonds

use ffdev_constants
use ffdev_variables
use ffdev_errors_dat

!===============================================================================

type, extends(ErrorFceType) :: TypeEFTBonds
    contains
        ! executive methods
        procedure   :: set_title_errfce         => ffdev_err_bonds_set_title
        procedure   :: calc_errfce              => ffdev_err_bonds_error
        procedure   :: print_individual_summary_errfce => ffdev_err_bonds_summary_targetset
        procedure   :: print_pts_summary_errfce => ffdev_err_bonds_summary
end type TypeEFTBonds

contains

!===============================================================================
! Subroutine:  ffdev_err_bonds_set_title
!===============================================================================

subroutine ffdev_err_bonds_set_title(err_item)

    implicit none
    class(TypeEFTBonds)    :: err_item
    ! --------------------------------------------------------------------------

    err_item%Title = 'Bonds'

end subroutine ffdev_err_bonds_set_title

! ==============================================================================
! subroutine ffdev_err_bonds_error
! opterr - calculate minimum for error optimization
! ==============================================================================

subroutine ffdev_err_bonds_error(err_item,opterr)

    use ffdev_targetset
    use ffdev_targetset_dat
    use ffdev_utils
    use ffdev_geometry
    use ffdev_errors_dat

    implicit none
    class(TypeEFTBonds) :: err_item
    logical             :: opterr
    ! --------------------------------------------
    integer             :: i,j,q,ai,aj
    real(DEVDP)         :: err,seterrbonds,totw
    real(DEVDP)         :: d0,dt
    ! --------------------------------------------------------------------------

    err_item%RepValue = 0.0d0
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

            do j=1,sets(i)%ngeos
                if( .not. sets(i)%geo(j)%trg_crd_optimized ) cycle

                d0 = ffdev_geometry_get_length(sets(i)%geo(j)%crd,ai,aj)
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

end subroutine ffdev_err_bonds_error

! ==============================================================================
! subroutine ffdev_err_bonds_summary_targetset
! ==============================================================================

subroutine ffdev_err_bonds_summary_targetset(err_item)

    use ffdev_geometry_utils

    implicit none
    class(TypeEFTBonds) :: err_item
    ! --------------------------------------------

    if( .not. err_item%PrintSummary ) return

    call ffdev_geometry_utils_targetset_stat_bonds(err_item%OnlyFFOpt)

end subroutine ffdev_err_bonds_summary_targetset

! ==============================================================================
! subroutine ffdev_err_bonds_summary
! ==============================================================================

subroutine ffdev_err_bonds_summary(err_item,top,geo,printsum)

    use ffdev_topology
    use ffdev_geometry
    use ffdev_geometry_utils

    implicit none
    class(TypeEFTBonds) :: err_item
    type(TOPOLOGY)      :: top
    type(GEOMETRY)      :: geo
    logical             :: printsum
    ! --------------------------------------------------------------------------

    if( .not. err_item%PrintSummary ) return
    if( .not. geo%trg_crd_optimized ) return

    if( printsum .eqv. .false. ) then
        printsum = top%nbonds .gt. 0
        return
    end if

    call ffdev_geometry_utils_comp_bonds(.false.,top,geo%trg_crd,geo%crd,err_item%OnlyFFOpt)

end subroutine ffdev_err_bonds_summary

! ------------------------------------------------------------------------------

end module ffdev_err_bonds
