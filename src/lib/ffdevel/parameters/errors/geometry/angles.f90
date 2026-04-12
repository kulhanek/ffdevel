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

module ffdev_err_angles

use ffdev_constants
use ffdev_variables
use ffdev_errors_dat

!===============================================================================

type, extends(ErrorFceType) :: TypeEFTAngles
    contains
        ! executive methods
        procedure   :: set_title_errfce         => ffdev_err_angles_set_title
        procedure   :: calc_errfce              => ffdev_err_angles_error
        procedure   :: print_individual_summary_errfce => ffdev_err_angles_summary_targetset
        procedure   :: print_pts_summary_errfce => ffdev_err_angles_summary
end type TypeEFTAngles

contains

!===============================================================================
! Subroutine:  ffdev_err_angles_set_title
!===============================================================================

subroutine ffdev_err_angles_set_title(err_item)

    implicit none
    class(TypeEFTAngles)    :: err_item
    ! --------------------------------------------------------------------------

    err_item%Title = 'Angles'

end subroutine ffdev_err_angles_set_title

! ==============================================================================
! subroutine ffdev_err_angles_error
! ==============================================================================

subroutine ffdev_err_angles_error(err_item,opterr)

    use ffdev_targetset_dat
    use ffdev_utils
    use ffdev_geometry
    use ffdev_errors_dat

    implicit none
    class(TypeEFTAngles)    :: err_item
    logical                 :: opterr
    ! --------------------------------------------
    integer                 :: i,j,q,ai,aj,ak
    real(DEVDP)             :: err,seterrangles,totw
    real(DEVDP)             :: d0,dt
    ! --------------------------------------------------------------------------

    err_item%RepValue = 0.0d0
    err_item%OptValue = 0.0d0
    if( .not. err_item%Enabled ) then
        if( opterr ) return
    end if

    seterrangles = 0.0
    totw = 0

    do i=1,nsets
        do q=1,sets(i)%top%nangles
            if( err_item%OnlyFFOpt ) then
                if( .not. sets(i)%top%angle_types(sets(i)%top%angles(q)%at)%ffoptactive ) cycle
            end if
            ai = sets(i)%top%angles(q)%ai
            aj = sets(i)%top%angles(q)%aj
            ak = sets(i)%top%angles(q)%ak

            do j=1,sets(i)%ngeos
                if( .not. sets(i)%geo(j)%trg_crd_optimized ) cycle

                d0 = ffdev_geometry_get_angle(sets(i)%geo(j)%crd,ai,aj,ak) * DEV_R2D
                dt = ffdev_geometry_get_angle(sets(i)%geo(j)%trg_crd,ai,aj,ak) * DEV_R2D
                err = d0 - dt
                seterrangles = seterrangles + sets(i)%geo(j)%weight * err**2
                totw = totw + sets(i)%geo(j)%weight
            end do
        end do
    end do

    if( totw .gt. 0 ) then
        err_item%RepValue = sqrt(seterrangles/totw)
    end if

    if( err_item%ScaleFac .gt. 0 ) then
        err_item%OptValue = err_item%Weight * err_item%RepValue**2 / err_item%ScaleFac**2
    end if

end subroutine ffdev_err_angles_error

! ==============================================================================
! subroutine ffdev_err_angles_summary_targetset
! ==============================================================================

subroutine ffdev_err_angles_summary_targetset(err_item)

    use ffdev_geometry_utils

    implicit none
    class(TypeEFTAngles) :: err_item
    ! --------------------------------------------

    if( .not. err_item%PrintSummary ) return

    call ffdev_geometry_utils_targetset_stat_angles(err_item%OnlyFFOpt)

end subroutine ffdev_err_angles_summary_targetset

! ==============================================================================
! subroutine ffdev_err_angles_summary
! ==============================================================================

subroutine ffdev_err_angles_summary(err_item,top,geo,printsum)

    use ffdev_topology
    use ffdev_geometry
    use ffdev_geometry_utils

    implicit none
    class(TypeEFTAngles)    :: err_item
    type(TOPOLOGY)          :: top
    type(GEOMETRY)          :: geo
    logical                 :: printsum
    ! --------------------------------------------------------------------------

    if( .not. err_item%PrintSummary ) return
    if( .not. geo%trg_crd_optimized ) return

    if( printsum .eqv. .false. ) then
        printsum = top%nangles .gt. 0
        return
    end if

    call ffdev_geometry_utils_comp_angles(.false.,top,geo%trg_crd,geo%crd,err_item%OnlyFFOpt)

end subroutine ffdev_err_angles_summary

! ------------------------------------------------------------------------------

end module ffdev_err_angles
