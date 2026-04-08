! ==============================================================================
! This file is part of FFDevel.
!    Copyright (C) 2026 Petr Kulhanek, kulhanek@chemi.muni.cz
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

module ffdev_err_ihess

use ffdev_constants
use ffdev_variables

contains

! ==============================================================================
! subroutine ffdev_err_ihess_init
! ==============================================================================

subroutine ffdev_err_ihess_init

    use ffdev_err_ihess_dat

    implicit none
    ! --------------------------------------------------------------------------

    EnableIHessError            = .false.
    PrintIHessErrorSummary      = .false.
    IHessErrorsWeightBonds      = 1.0
    IHessErrorsWeightAngles     = 1.0
    IHessErrorsExcludeDihedrals = .true.
    IHessErrorsExcludeImpropers = .true.
    IHessErrorsExcludeNB        = .true.
    OnlyFFOptIHess              = .false.

end subroutine ffdev_err_ihess_init

! ==============================================================================
! subroutine ffdev_err_ihess_error
! ==============================================================================

subroutine ffdev_err_ihess_error(error)

    use ffdev_targetset
    use ffdev_targetset_dat
    use ffdev_utils
    use ffdev_geometry
    use ffdev_errors_dat
    use ffdev_err_ihess_dat

    implicit none
    type(FFERROR_TYPE)  :: error
    ! --------------------------------------------
    integer             :: s,i,j,bt,at,nihess_b,nihess_a
    real(DEVDP)         :: k1,k2,ihess_bse,ihess_ase,diff
    ! --------------------------------------------------------------------------

    error%ihess_bonds  = 0.0d0
    error%ihess_angles = 0.0d0

    ! calculate ihess
    do s=1,nsets
        do i=1,sets(s)%ngeos
            call ffdev_err_ihess_error_geo(sets(s)%top,sets(s)%geo(i))
        end do

    end do

    ! calculate error
    nihess_b = 0
    ihess_bse = 0.0

    do s=1,nsets
        do i=1,sets(s)%top%nbonds
            if( OnlyFFOptIHess) then
                if( .not. sets(s)%top%bond_types(sets(s)%top%bonds(i)%bt)%ffoptactive ) cycle
            end if
            bt = sets(s)%top%bonds(i)%bt
            k1 = sets(s)%top%bond_types(bt)%k
            do j=1,sets(s)%ngeos
                k2 = sets(s)%geo(j)%trg_ihess_bonds(i)
                diff = k2 - k1
                ihess_bse = ihess_bse + diff**2
                nihess_b = nihess_b + 1
            end do
        end do
    end do

    if( nihess_b .gt. 0 ) then
        error%ihess_bonds = sqrt(ihess_bse / real(nihess_b))
    end if

    nihess_a = 0
    ihess_ase = 0.0

    do s=1,nsets
        do i=1,sets(s)%top%nangles
            if( OnlyFFOptIHess) then
                if( .not. sets(s)%top%angle_types(sets(s)%top%angles(i)%at)%ffoptactive ) cycle
            end if
            at = sets(s)%top%angles(i)%at
            k1 = sets(s)%top%angle_types(at)%k
            do j=1,sets(s)%ngeos
                k2 = sets(s)%geo(j)%trg_ihess_angles(i)
                diff = k2 - k1
                ihess_ase = ihess_ase + diff**2
                nihess_a = nihess_a + 1
            end do
        end do
    end do

    if( nihess_a .gt. 0 ) then
        error%ihess_angles = sqrt(ihess_ase / real(nihess_a))
    end if

end subroutine ffdev_err_ihess_error

! ==============================================================================
! subroutine ffdev_err_ihess_error_geo
! ==============================================================================

subroutine ffdev_err_ihess_error_geo(top,geo)

    use ffdev_topology
    use ffdev_geometry
    use ffdev_geometry_utils
    use ffdev_gradient_utils
    use ffdev_hessian
    use ffdev_hessian_utils
    use ffdev_utils
    use ffdev_err_ihess_dat
    use ffdev_nbmode_LJ

    implicit none
    type(TOPOLOGY)  :: top
    type(GEOMETRY)  :: geo
    ! --------------------------------------------
    type(GEOMETRY)  :: geo_backup
    ! -----------------------------------------------------------------------------

    if( .not. associated(geo%trg_ihess) ) then
        call ffdev_hessian_allocate_trg_ihess(top,geo)
    end if

    if( IHessErrorsExcludeDihedrals .or. IHessErrorsExcludeImpropers .or. IHessErrorsExcludeNB ) then
        call ffdev_geometry_init(geo_backup)
        call ffdev_geometry_allocate(geo_backup,geo%natoms)
        ! copy target coordinates
        geo_backup%crd = geo%trg_crd
        call ffdev_gradient_allocate(geo_backup)
        geo_backup%grd = 0.0d0
        call ffdev_hessian_allocate(geo_backup)
        geo_backup%hess = 0.0d0
    end if

    if( IHessErrorsExcludeDihedrals ) then
        call ffdev_hessian_dihedrals(top,geo_backup)
    end if

    if( IHessErrorsExcludeImpropers ) then
        call ffdev_hessian_impropers(top,geo_backup)
    end if

    if( IHessErrorsExcludeNB ) then
        call ffdev_hessian_nb_lj(top,geo_backup)
    end if

    if( IHessErrorsExcludeDihedrals .or. IHessErrorsExcludeImpropers .or. IHessErrorsExcludeNB ) then
        geo%trg_ihess = geo%trg_hess - geo_backup%hess
        call ffdev_geometry_destroy(geo_backup)
    else
        geo%trg_ihess = geo%trg_hess
    end if

    ! calculate ihess
    call ffdev_hessian_calc_trg_ihess(top,geo)

end subroutine ffdev_err_ihess_error_geo

! ==============================================================================
! subroutine ffdev_err_ihess_summary
! ==============================================================================

subroutine ffdev_err_ihess_summary(top,geo,printsum)

    use ffdev_topology
    use ffdev_geometry
    use ffdev_geometry_utils
    use ffdev_hessian_utils

    implicit none
    type(TOPOLOGY)  :: top
    type(GEOMETRY)  :: geo
    logical         :: printsum
    ! --------------------------------------------------------------------------

    if( printsum .eqv. .false. ) then
        printsum = (top%nbonds .gt. 0) .or. (top%nangles .le. 0)
        return
    end if

    ! print results
    call ffdev_hessian_print_trg_ihess_bonds(top,geo)
    call ffdev_hessian_print_trg_ihess_angles(top,geo)

end subroutine ffdev_err_ihess_summary

! ------------------------------------------------------------------------------

end module ffdev_err_ihess
