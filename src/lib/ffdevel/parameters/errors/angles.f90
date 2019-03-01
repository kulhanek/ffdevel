! ==============================================================================
! This file is part of FFDevel.
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

contains

! ==============================================================================
! subroutine ffdev_err_angles_init
! ==============================================================================

subroutine ffdev_err_angles_init

    use ffdev_err_angles_dat

    implicit none
    ! --------------------------------------------------------------------------

    EnableAngleError         = .false.
    PrintAngleErrorSummary   = .false.
    AngleErrorWeight         = DEV_D2R

end subroutine ffdev_err_angles_init

! ==============================================================================
! subroutine ffdev_err_angles_error
! ==============================================================================

subroutine ffdev_err_angles_error(error)

    use ffdev_targetset_dat
    use ffdev_utils   
    use ffdev_geometry
    use ffdev_errors_dat

    implicit none
    type(FFERROR_TYPE)  :: error
    ! --------------------------------------------
    integer             :: i,j,q,nangles,ai,aj,ak
    real(DEVDP)         :: err,seterrangles
    real(DEVDP)         :: d0,dt
    ! --------------------------------------------------------------------------

    error%angles = 0.0

    seterrangles = 0.0
    nangles = 0

    do i=1,nsets
        do j=1,sets(i)%ngeos
            ! ------------------------------------------------------------------
            if( sets(i)%geo(j)%trg_crd_optimized ) then
                do q=1,sets(i)%top%nangles
                    ai = sets(i)%top%angles(q)%ai
                    aj = sets(i)%top%angles(q)%aj
                    ak = sets(i)%top%angles(q)%ak
                    d0 = ffdev_geometry_get_angle(sets(i)%geo(j)%crd,ai,aj,ak) * DEV_R2D
                    dt = ffdev_geometry_get_angle(sets(i)%geo(j)%trg_crd,ai,aj,ak) * DEV_R2D
                    nangles = nangles + 1
                    err = d0 - dt
                    seterrangles = seterrangles + sets(i)%geo(j)%weight * err**2
                end do
            end if          
        end do
    end do

    if( nangles .gt. 0 ) then
        error%angles = sqrt(seterrangles/real(nangles))
    end if 

end subroutine ffdev_err_angles_error

! ==============================================================================
! subroutine ffdev_err_angles_summary
! ==============================================================================

subroutine ffdev_err_angles_summary(top,geo)

    use ffdev_topology
    use ffdev_geometry
    use ffdev_geometry_utils

    implicit none
    type(TOPOLOGY)     :: top
    type(GEOMETRY)     :: geo
    ! --------------------------------------------------------------------------

    call ffdev_geometry_utils_comp_angles(.false.,top,geo%trg_crd,geo%crd)

end subroutine ffdev_err_angles_summary

! ------------------------------------------------------------------------------

end module ffdev_err_angles
