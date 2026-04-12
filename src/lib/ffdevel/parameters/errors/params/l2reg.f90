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

module ffdev_err_l2reg

use ffdev_constants
use ffdev_variables
use ffdev_errors_dat

!===============================================================================

type, extends(ErrorFceType) :: TypeEFTL2Reg

    integer         :: Realm     ! parameter realm

    contains
        ! executive methods
        procedure   :: init_errfce              => ffdev_err_l2reg_init
        procedure   :: load_errfce              => ffdev_err_l2reg_ctrl
        procedure   :: set_title_errfce         => ffdev_err_l2reg_set_title
        procedure   :: calc_errfce              => ffdev_err_l2reg_error
end type TypeEFTL2Reg

contains

! ==============================================================================
! subroutine ffdev_err_l2reg_init
! ==============================================================================

subroutine ffdev_err_l2reg_init(err_item)

    use ffdev_parameters_dat

    implicit none
    class(TypeEFTL2Reg)    :: err_item
    ! --------------------------------------------------------------------------

    call err_item%ErrorFceType%init_errfce()

    err_item%Realm          = REALM_DIH_C
    err_item%ScaleFac       = 0.0d0 ! auto

end subroutine ffdev_err_l2reg_init

! ==============================================================================
! subroutine ffdev_err_l2reg_ctrl
! ==============================================================================

subroutine ffdev_err_l2reg_ctrl(err_item,fin)

    use ffdev_utils
    use prmfile
    use ffdev_parameters_dat
    use ffdev_parameters

    implicit none
    class(TypeEFTL2Reg)     :: err_item
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

160  format ('Parameter realm (realm)                = ',a12)
165  format ('Parameter realm (realm)                = ',a12,'                  (default)')

end subroutine ffdev_err_l2reg_ctrl

!===============================================================================
! Subroutine:  ffdev_err_l2reg_set_title
!===============================================================================

subroutine ffdev_err_l2reg_set_title(err_item)

    use ffdev_parameters
    use ffdev_utils

    implicit none
    class(TypeEFTL2Reg)    :: err_item
    ! --------------------------------------------------------------------------

    err_item%Title = 'L2(' // trim(ffdev_parameters_get_realm_name(err_item%Realm)) // ')'

end subroutine ffdev_err_l2reg_set_title

! ==============================================================================
! subroutine ffdev_err_l2reg_error
! ==============================================================================

subroutine ffdev_err_l2reg_error(err_item,opterr)

    use ffdev_targetset
    use ffdev_targetset_dat
    use ffdev_utils
    use ffdev_geometry
    use ffdev_errors_dat
    use ffdev_parameters_dat
    use ffdev_parameters

    implicit none
    class(TypeEFTL2Reg) :: err_item
    logical             :: opterr
    ! --------------------------------------------
    integer             :: i, pc
    real(DEVDP)         :: v, used_scale_fac
    ! --------------------------------------------------------------------------

    err_item%RepValue = 0.0d0
    err_item%OptValue = 0.0d0
    if( .not. err_item%Enabled ) then
        if( opterr ) return
    end if

    pc = 0
    do i=1,nparams
        ! skip different realms
        if( params(i)%realm .ne. err_item%Realm ) cycle

        ! check if it is activated
        if( err_item%OnlyFFOpt ) then
            if( .not. params(i)%enabled ) cycle
        end if

        ! calculate error
        v = params(i)%value
        err_item%RepValue = err_item%RepValue + v**2
        pc = pc + 1
    end do

    used_scale_fac = err_item%ScaleFac
    if( used_scale_fac .eq. 0 ) then
        used_scale_fac = pc * max( abs(ffdev_params_get_lower_bound(err_item%Realm)), &
                              abs(ffdev_params_get_upper_bound(err_item%Realm)) )**2
    end if

    if( used_scale_fac .gt. 0 ) then
        err_item%OptValue = err_item%Weight * err_item%RepValue / used_scale_fac
    end if

end subroutine ffdev_err_l2reg_error

! ------------------------------------------------------------------------------

end module ffdev_err_l2reg
