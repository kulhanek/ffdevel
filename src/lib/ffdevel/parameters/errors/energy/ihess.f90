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
use ffdev_errors_dat

!===============================================================================

type, extends(ErrorFceType) :: TypeEFTIHess

    integer         :: Realm                ! BOND_K or ANGLE_K
    logical         :: ExcludeDihedrals
    logical         :: ExcludeImpropers
    logical         :: ExcludeNB

    contains
        ! executive methods
        procedure   :: init_errfce              => ffdev_err_ihess_init
        procedure   :: load_errfce              => ffdev_err_ihess_ctrl
        procedure   :: set_title_errfce         => ffdev_err_ihess_set_title
        procedure   :: calc_errfce              => ffdev_err_ihess_error
        procedure   :: print_pts_summary_errfce => ffdev_err_ihess_summary
end type TypeEFTIHess

contains

! ==============================================================================
! subroutine ffdev_err_ihess_init
! ==============================================================================

subroutine ffdev_err_ihess_init(err_item)

    use ffdev_parameters_dat

    implicit none
    class(TypeEFTIHess)    :: err_item
    ! --------------------------------------------------------------------------

    call err_item%ErrorFceType%init_errfce()

    err_item%Realm            = REALM_BOND_K
    err_item%ExcludeDihedrals = .true.
    err_item%ExcludeImpropers = .true.
    err_item%ExcludeNB        = .true.

end subroutine ffdev_err_ihess_init

! ==============================================================================
! subroutine ffdev_err_ihess_ctrl
! ==============================================================================

subroutine ffdev_err_ihess_ctrl(err_item,fin)

    use ffdev_utils
    use prmfile
    use ffdev_parameters_dat
    use ffdev_parameters

    implicit none
    class(TypeEFTIHess)     :: err_item
    type(PRMFILE_TYPE)      :: fin
    ! --------------------------------------------
    character(PRMFILE_MAX_PATH) :: realm
    ! --------------------------------------------------------------------------

    call err_item%ErrorFceType%load_errfce(fin)

    realm = ffdev_parameters_get_realm_name(err_item%Realm)

    if( prmfile_get_string_by_key(fin,'realm', realm)) then
        write(DEV_OUT,160) realm
    else
        write(DEV_OUT,165) realm
    end if

    err_item%Realm = ffdev_parameters_get_realmid(realm)

    if( prmfile_get_logical_by_key(fin,'exclude_dihedrals', err_item%ExcludeDihedrals)) then
        write(DEV_OUT,170) prmfile_onoff(err_item%ExcludeDihedrals)
    else
        write(DEV_OUT,175) prmfile_onoff(err_item%ExcludeDihedrals)
    end if
    if( prmfile_get_logical_by_key(fin,'exclude_impropers', err_item%ExcludeImpropers)) then
        write(DEV_OUT,180) prmfile_onoff(err_item%ExcludeImpropers)
    else
        write(DEV_OUT,185) prmfile_onoff(err_item%ExcludeImpropers)
    end if
    if( prmfile_get_logical_by_key(fin,'exclude_nb', err_item%ExcludeNB)) then
        write(DEV_OUT,190) prmfile_onoff(err_item%ExcludeNB)
    else
        write(DEV_OUT,195) prmfile_onoff(err_item%ExcludeNB)
    end if

160  format ('Realm (realm)                          = ',a12)
165  format ('Realm (realm)                          = ',a12,'                  (default)')

170  format ('Exclude dihedrals (exlude_dihedrals)   = ',a12)
175  format ('Exclude dihedrals (exlude_dihedrals)   = ',a12,'                  (default)')

180  format ('Exclude impropers (exlude_impropers)   = ',a12)
185  format ('Exclude impropers (exlude_impropers)   = ',a12,'                  (default)')

190  format ('Exclude NB (exlude_nb)                 = ',a12)
195  format ('Exclude NB (exlude_nb)                 = ',a12,'                  (default)')

end subroutine ffdev_err_ihess_ctrl

!===============================================================================
! Subroutine:  ffdev_err_ihess_set_title
!===============================================================================

subroutine ffdev_err_ihess_set_title(err_item)

    use ffdev_parameters_dat
    use ffdev_utils

    implicit none
    class(TypeEFTIHess)    :: err_item
    ! --------------------------------------------------------------------------

    select case(err_item%Realm)
        case(REALM_BOND_K)
            err_item%Title = 'IH(Bonds)'
        case(REALM_ANGLE_K)
            err_item%Title = 'IH(Angles)'
        case default
            call ffdev_utils_exit(DEV_ERR,1,'Unsupported realm in ffdev_err_ihess_set_title!')
    end select

end subroutine ffdev_err_ihess_set_title

! ==============================================================================
! subroutine ffdev_err_ihess_error
! ==============================================================================

subroutine ffdev_err_ihess_error(err_item,opterr)

    use ffdev_targetset
    use ffdev_targetset_dat
    use ffdev_utils
    use ffdev_geometry
    use ffdev_errors_dat
    use ffdev_parameters_dat

    implicit none
    class(TypeEFTIHess)     :: err_item
    logical                 :: opterr
    ! --------------------------------------------
    integer             :: s,i,j,bt,at,nihess_b,nihess_a
    real(DEVDP)         :: k1,k2,ihess_bse,ihess_ase,diff
    ! --------------------------------------------------------------------------

    err_item%RepValue = 0.0d0
    if( .not. err_item%Enabled ) then
        if( opterr ) return
    end if

    ! calculate ihess
    do s=1,nsets
        do i=1,sets(s)%ngeos
            call ffdev_err_ihess_error_geo(err_item,sets(s)%top,sets(s)%geo(i))
        end do

    end do

    ! calculate error
    select case(err_item%Realm)
        case(REALM_BOND_K)
            nihess_b = 0
            ihess_bse = 0.0

            do s=1,nsets
                do i=1,sets(s)%top%nbonds
                    if( err_item%OnlyFFOpt ) then
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
                err_item%RepValue = sqrt(ihess_bse / real(nihess_b))
            end if
        case(REALM_ANGLE_K)
            nihess_a = 0
            ihess_ase = 0.0

            do s=1,nsets
                do i=1,sets(s)%top%nangles
                    if( err_item%OnlyFFOpt ) then
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
                err_item%RepValue = sqrt(ihess_ase / real(nihess_a))
            end if
        case default
            call ffdev_utils_exit(DEV_ERR,1,'Unsupported realm in ffdev_err_ihess_error!')
    end select

end subroutine ffdev_err_ihess_error

! ==============================================================================
! subroutine ffdev_err_ihess_error_geo
! ==============================================================================

subroutine ffdev_err_ihess_error_geo(err_item,top,geo)

    use ffdev_topology
    use ffdev_geometry
    use ffdev_geometry_utils
    use ffdev_gradient_utils
    use ffdev_hessian
    use ffdev_hessian_utils
    use ffdev_utils
    use ffdev_nbmode_LJ

    implicit none
    class(TypeEFTIHess) :: err_item
    type(TOPOLOGY)      :: top
    type(GEOMETRY)      :: geo
    ! --------------------------------------------
    type(GEOMETRY)      :: geo_backup
    ! -----------------------------------------------------------------------------

    if( .not. associated(geo%trg_ihess) ) then
        call ffdev_hessian_allocate_trg_ihess(top,geo)
    end if

    if( err_item%ExcludeDihedrals .or. err_item%ExcludeImpropers .or. err_item%ExcludeNB ) then
        call ffdev_geometry_init(geo_backup)
        call ffdev_geometry_allocate(geo_backup,geo%natoms)
        ! copy target coordinates
        geo_backup%crd = geo%trg_crd
        call ffdev_gradient_allocate(geo_backup)
        geo_backup%grd = 0.0d0
        call ffdev_hessian_allocate(geo_backup)
        geo_backup%hess = 0.0d0
    end if

    if( err_item%ExcludeDihedrals ) then
        call ffdev_hessian_dihedrals(top,geo_backup)
    end if

    if( err_item%ExcludeImpropers ) then
        call ffdev_hessian_impropers(top,geo_backup)
    end if

    if( err_item%ExcludeNB ) then
        call ffdev_hessian_nb_lj(top,geo_backup)
    end if

    if( err_item%ExcludeDihedrals .or. err_item%ExcludeImpropers .or. err_item%ExcludeNB ) then
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

subroutine ffdev_err_ihess_summary(err_item,top,geo,printsum)

    use ffdev_topology
    use ffdev_geometry
    use ffdev_geometry_utils
    use ffdev_hessian_utils
    use ffdev_parameters_dat
    use ffdev_utils

    implicit none
    class(TypeEFTIHess) :: err_item
    type(TOPOLOGY)      :: top
    type(GEOMETRY)      :: geo
    logical             :: printsum
    ! --------------------------------------------------------------------------

    if( .not. err_item%PrintSummary ) return

    if( printsum .eqv. .false. ) then
        printsum = (top%nbonds .gt. 0) .or. (top%nangles .le. 0)
        return
    end if

    ! print results
    select case(err_item%Realm)
        case(REALM_BOND_K)
            call ffdev_hessian_print_trg_ihess_bonds(top,geo,err_item%OnlyFFOpt)
        case(REALM_ANGLE_K)
            call ffdev_hessian_print_trg_ihess_angles(top,geo,err_item%OnlyFFOpt)
        case default
            call ffdev_utils_exit(DEV_ERR,1,'Unsupported realm in ffdev_err_ihess_summary!')
    end select

end subroutine ffdev_err_ihess_summary

! ------------------------------------------------------------------------------

end module ffdev_err_ihess
