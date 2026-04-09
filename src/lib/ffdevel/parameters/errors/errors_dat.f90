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

module ffdev_errors_dat

use ffdev_constants
use ffdev_variables

! ------------------------------------------------------------------------------

integer,parameter       :: EE_ABS       = 1 ! absolute
integer,parameter       :: EE_REL       = 2 ! relative
integer,parameter       :: EE_LOG       = 3 ! log


! base class for error function ------------------------------------------------
type ErrorFceType
    ! setup
    logical                     :: Enabled      ! enable the error as the part of the optimized objective function
    logical                     :: PrintSummary ! print error summary
    logical                     :: OnlyFFOpt    ! consider items that are related to optimized parameters
    real(DEVDP)                 :: Weight       ! weight of error function

    ! calculated
    character(MAX_TITLE)        :: Title        ! short description

    ! value
    real(DEVDP)                 :: ErrFceValue  ! value of error function

    contains
        ! executive methods
        procedure               :: init_errfce
        procedure               :: load_errfce
        procedure               :: set_title_errfce
        procedure               :: setup_domains_errfce
        procedure               :: calc_errfce
        procedure               :: print_individual_summary_errfce
        procedure               :: print_set_summary_errfce
        procedure               :: print_pts_summary_errfce
end type ErrorFceType

! ------------------------------------------------------------------------------

type ErrorFcePointer
    class(ErrorFceType),pointer   :: ErrFce
end type ErrorFcePointer

! ------------------------------------------------------------------------------

integer                             :: NumOfErrorFces = 0   ! number of error functions
type(ErrorFcePointer),allocatable   :: ErrorFceList(:)      ! list of error functions

! ------------------------------------------------------------------------------

real(DEVDP)             :: FFError

! ------------------------------------------------------------------------------

integer, parameter      :: SMMLOG_INITIAL       = 1
integer, parameter      :: SMMLOG_INTERMEDIATE  = 2
integer, parameter      :: SMMLOG_FINAL         = 3
integer, parameter      :: SMMLOG_BEST          = 4

! ------------------------------------------------------------------------------
! internal setup
logical                 :: errors_calc_ene      = .false.
logical                 :: errors_calc_sapt     = .false.
logical                 :: errors_calc_grad     = .false.
logical                 :: errors_calc_hess     = .false.

contains
!===============================================================================
! Subroutine:  init_errfce
!===============================================================================

subroutine init_errfce(err_item)

    implicit none
    class(ErrorFceType)    :: err_item
    ! --------------------------------------------------------------------------

    err_item%Title           = ''

    err_item%Enabled         = .false.
    err_item%PrintSummary    = .false.
    err_item%OnlyFFOpt       = .false.
    err_item%Weight          = 1.0d0

    err_item%ErrFceValue     = 0.0d0

end subroutine init_errfce

!===============================================================================
! Subroutine:  load_errfce
!===============================================================================

subroutine load_errfce(err_item,fin)

    use prmfile

    implicit none
    class(ErrorFceType) :: err_item
    type(PRMFILE_TYPE)  :: fin
    ! --------------------------------------------------------------------------

    ! load setup
    if( prmfile_get_logical_by_key(fin,'enabled', err_item%Enabled)) then
        write(DEV_OUT,110) prmfile_onoff(err_item%Enabled)
    else
        write(DEV_OUT,115) prmfile_onoff(err_item%Enabled)
    end if
    if( prmfile_get_logical_by_key(fin,'summary', err_item%PrintSummary)) then
        write(DEV_OUT,130) prmfile_onoff( err_item%PrintSummary)
    else
        write(DEV_OUT,135) prmfile_onoff( err_item%PrintSummary)
    end if
    if( prmfile_get_real8_by_key(fin,'weight', err_item%Weight)) then
        write(DEV_OUT,120) err_item%Weight
    else
        write(DEV_OUT,125) err_item%Weight
    end if
    if( prmfile_get_logical_by_key(fin,'onlyffopt', err_item%OnlyFFOpt)) then
        write(DEV_OUT,140) prmfile_onoff(err_item%OnlyFFOpt)
    else
        write(DEV_OUT,145) prmfile_onoff(err_item%OnlyFFOpt)
    end if

110  format ('Error enabled (enabled)                = ',a12)
115  format ('Error enabled (enabled)                = ',a12,'                  (default)')
130  format ('Print error summary (summary)          = ',a12)
135  format ('Print error summary (summary)          = ',a12,'                  (default)')
120  format ('Error weight (weight)                  = ',f21.8)
125  format ('Error weight (weight)                  = ',f21.8,'         (default)')
140  format ('Only FFopt related (onlyffopt)         = ',a12)
145  format ('Only FFopt related (onlyffopt)         = ',a12,'                  (default)')

end subroutine load_errfce

!===============================================================================
! Subroutine:  set_title_errfce
!===============================================================================

subroutine set_title_errfce(err_item)

    implicit none
    class(ErrorFceType)    :: err_item
    ! --------------------------------------------------------------------------

    ! disable unused variable warning
    ignored_arg__ = same_type_as(err_item,err_item)

end subroutine set_title_errfce

!===============================================================================
! Subroutine:  setup_domains_errfce
!===============================================================================

subroutine setup_domains_errfce(err_item,opterr)

    implicit none
    class(ErrorFceType) :: err_item
    logical             :: opterr
    ! --------------------------------------------------------------------------

    ! disable unused variable warning
    ignored_arg__ = same_type_as(err_item,err_item)
    ignored_arg__ = opterr .eqv. opterr

end subroutine setup_domains_errfce

!===============================================================================
! Subroutine:  calc_errfce
! opterr - calculate minimum for error optimization
!===============================================================================

subroutine calc_errfce(err_item,opterr)

    implicit none
    class(ErrorFceType) :: err_item
    logical             :: opterr
    ! --------------------------------------------------------------------------

    err_item%ErrFceValue = 0.0d0

    ! disable unused variable warning
    ignored_arg__ = opterr .eqv. opterr

end subroutine calc_errfce

!===============================================================================
! Subroutine:  print_individual_summary_errfce
!===============================================================================

subroutine print_individual_summary_errfce(err_item)

    implicit none
    class(ErrorFceType) :: err_item
    ! --------------------------------------------------------------------------

    ! disable unused variable warning
    ignored_arg__ = same_type_as(err_item,err_item)

end subroutine print_individual_summary_errfce

!===============================================================================
! Subroutine:  print_set_summary_errfce
!===============================================================================

subroutine print_set_summary_errfce(err_item,set,printsum)

    use ffdev_targetset_dat

    implicit none
    class(ErrorFceType) :: err_item
    type(TARGETSET)     :: set
    logical             :: printsum
    ! --------------------------------------------------------------------------

    ! disable unused variable warning
    ignored_arg__ = same_type_as(err_item,err_item)
    ignored_arg__ = same_type_as(set,set)
    ignored_arg__ = printsum .eqv. printsum

end subroutine print_set_summary_errfce

!===============================================================================
! Subroutine:  print_pts_summary_errfce
!===============================================================================

subroutine print_pts_summary_errfce(err_item,top,geo,printsum)

    use ffdev_topology_dat
    use ffdev_geometry_dat

    implicit none
    class(ErrorFceType) :: err_item
    type(TOPOLOGY)      :: top
    type(GEOMETRY)      :: geo
    logical             :: printsum
    ! --------------------------------------------------------------------------

    ! disable unused variable warning
    ignored_arg__ = same_type_as(err_item,err_item)
    ignored_arg__ = same_type_as(top,top)
    ignored_arg__ = same_type_as(geo,geo)
    ignored_arg__ = printsum .eqv. printsum

end subroutine print_pts_summary_errfce

! ------------------------------------------------------------------------------

end module ffdev_errors_dat
