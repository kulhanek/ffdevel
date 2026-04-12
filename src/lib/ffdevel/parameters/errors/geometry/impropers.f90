! ==============================================================================
! This file is part of FFDevel.
!    Copyright (C) 2026 Petr Kulhanek, kulhanek@chemi.muni.cz
!    Copyright (C) 2019 Petr Kulhanek, kulhanek@chemi.muni.cz
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

module ffdev_err_impropers

use ffdev_constants
use ffdev_variables
use ffdev_errors_dat

!===============================================================================

type, extends(ErrorFceType) :: TypeEFTImpropers

    logical         :: LockToPhase

    contains
        ! executive methods
        procedure   :: init_errfce              => ffdev_err_impropers_init
        procedure   :: load_errfce              => ffdev_err_impropers_ctrl
        procedure   :: set_title_errfce         => ffdev_err_impropers_set_title
        procedure   :: calc_errfce              => ffdev_err_impropers_error
        procedure   :: print_pts_summary_errfce => ffdev_err_impropers_summary
end type TypeEFTImpropers

contains

!===============================================================================
! Subroutine:  init_errfce
!===============================================================================

subroutine ffdev_err_impropers_init(err_item)

    implicit none
    class(TypeEFTImpropers) :: err_item
    ! --------------------------------------------------------------------------

    call err_item%ErrorFceType%init_errfce()

    err_item%LockToPhase = .false.

end subroutine ffdev_err_impropers_init

! ==============================================================================
! subroutine ffdev_err_impropers_ctrl
! ==============================================================================

subroutine ffdev_err_impropers_ctrl(err_item,fin)

    use ffdev_utils
    use prmfile

    implicit none
    class(TypeEFTImpropers) :: err_item
    type(PRMFILE_TYPE)      :: fin
    ! --------------------------------------------------------------------------

    call err_item%ErrorFceType%load_errfce(fin)

    if( prmfile_get_logical_by_key(fin,'lock2phase', err_item%LockToPhase)) then
        write(DEV_OUT,10) prmfile_onoff(err_item%LockToPhase)
    else
        write(DEV_OUT,15) prmfile_onoff(err_item%LockToPhase)
    end if

10  format ('Lock to phase angle (lock2phase)       = ',a12)
15  format ('Lock to phase angle (lock2phase)       = ',a12,'                  (default)')

end subroutine ffdev_err_impropers_ctrl

!===============================================================================
! Subroutine:  ffdev_err_impropers_set_title
!===============================================================================

subroutine ffdev_err_impropers_set_title(err_item)

    implicit none
    class(TypeEFTImpropers)    :: err_item
    ! --------------------------------------------------------------------------

    err_item%Title = 'Impropers'

end subroutine ffdev_err_impropers_set_title

! ==============================================================================
! subroutine ffdev_err_impropers_error
! ==============================================================================

subroutine ffdev_err_impropers_error(err_item,opterr)

    use ffdev_targetset_dat
    use ffdev_utils
    use ffdev_geometry
    use ffdev_errors_dat

    implicit none
    class(TypeEFTImpropers) :: err_item
    logical                 :: opterr
    ! --------------------------------------------
    integer             :: i,j,q,ai,aj,ak,al,idt
    real(DEVDP)         :: err,seterrimpropers,totw
    real(DEVDP)         :: d0,dt
    ! --------------------------------------------------------------------------

    err_item%RepValue = 0.0d0
    if( .not. err_item%Enabled ) then
        if( opterr ) return
    end if

    seterrimpropers = 0.0
    totw = 0

    do i=1,nsets
        do q=1,sets(i)%top%nimpropers
            if( err_item%OnlyFFOpt ) then
                if( .not. sets(i)%top%improper_types(sets(i)%top%impropers(q)%dt)%ffoptactive ) cycle
            end if
            ai = sets(i)%top%impropers(q)%ai
            aj = sets(i)%top%impropers(q)%aj
            ak = sets(i)%top%impropers(q)%ak
            al = sets(i)%top%impropers(q)%al

            do j=1,sets(i)%ngeos
                if( .not. sets(i)%geo(j)%trg_crd_optimized ) cycle

                d0 = ffdev_geometry_get_improper(sets(i)%geo(j)%crd,ai,aj,ak,al)
                if( err_item%LockToPhase ) then
                    idt = sets(i)%top%impropers(q)%dt
                    dt = sets(i)%top%improper_types(idt)%g
                else
                    dt = ffdev_geometry_get_improper(sets(i)%geo(j)%trg_crd,ai,aj,ak,al)
                end if
                err = ffdev_geometry_get_dihedral_deviation(d0,dt) ! this needs values in RAD
                err = err * DEV_R2D
                seterrimpropers = seterrimpropers + sets(i)%geo(j)%weight * err**2
                totw = totw + sets(i)%geo(j)%weight
            end do
        end do
    end do

    if( totw .gt. 0 ) then
        err_item%RepValue = sqrt(seterrimpropers/totw)
    end if

end subroutine ffdev_err_impropers_error

! ==============================================================================
! subroutine ffdev_err_impropers_summary
! ==============================================================================

subroutine ffdev_err_impropers_summary(err_item,top,geo,printsum)

    use ffdev_topology
    use ffdev_geometry
    use ffdev_geometry_utils

    implicit none
    class(TypeEFTImpropers) :: err_item
    type(TOPOLOGY)          :: top
    type(GEOMETRY)          :: geo
    logical                 :: printsum
    ! --------------------------------------------------------------------------

    if( .not. err_item%PrintSummary ) return
    if( .not. geo%trg_crd_optimized ) return

    if( printsum .eqv. .false. ) then
        printsum = top%nimpropers .gt. 0
        return
    end if

    call ffdev_geometry_utils_comp_impropers(.false.,top,geo%trg_crd,geo%crd,err_item%LockToPhase)

end subroutine ffdev_err_impropers_summary

! ------------------------------------------------------------------------------

end module ffdev_err_impropers
