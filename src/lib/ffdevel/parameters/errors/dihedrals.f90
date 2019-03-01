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

module ffdev_err_dihedrals

use ffdev_constants

contains

! ==============================================================================
! subroutine ffdev_err_dihedrals_init
! ==============================================================================

subroutine ffdev_err_dihedrals_init

    use ffdev_err_dihedrals_dat

    implicit none
    ! --------------------------------------------------------------------------

    EnableDihedralError           = .false.
    PrintDihedralErrorSummary     = .false.
    DihedralErrorWeight           = DEV_D2R

end subroutine ffdev_err_dihedrals_init

! ==============================================================================
! subroutine ffdev_err_dihedrals_error
! ==============================================================================

subroutine ffdev_err_dihedrals_error(error)

    use ffdev_targetset_dat
    use ffdev_utils   
    use ffdev_geometry
    use ffdev_errors_dat

    implicit none
    type(FFERROR_TYPE)  :: error
    ! --------------------------------------------
    integer             :: i,j,q,ndihedrals,ai,aj,ak,al
    real(DEVDP)         :: err,seterrdihedrals
    real(DEVDP)         :: d0,dt
    ! --------------------------------------------------------------------------

    error%dihedrals = 0.0

    seterrdihedrals = 0.0
    ndihedrals = 0

    if( sets(i)%geo(j)%trg_crd_optimized ) then
        do q=1,sets(i)%top%ndihedrals
            ai = sets(i)%top%dihedrals(q)%ai
            aj = sets(i)%top%dihedrals(q)%aj
            ak = sets(i)%top%dihedrals(q)%ak
            al = sets(i)%top%dihedrals(q)%al
            d0 = ffdev_geometry_get_dihedral(sets(i)%geo(j)%crd,ai,aj,ak,al) * DEV_R2D
            dt = ffdev_geometry_get_dihedral(sets(i)%geo(j)%trg_crd,ai,aj,ak,al) * DEV_R2D
            ndihedrals = ndihedrals + 1
            err = ffdev_geometry_get_dihedral_deviation(d0,dt)
            seterrdihedrals = seterrdihedrals + sets(i)%geo(j)%weight * err**2
        end do
    end if

    if( ndihedrals .gt. 0 ) then
        error%dihedrals = sqrt(seterrdihedrals/real(ndihedrals))
    end if 

end subroutine ffdev_err_dihedrals_error

! ==============================================================================
! subroutine ffdev_err_dihedrals_summary
! ==============================================================================

subroutine ffdev_err_dihedrals_summary(top,geo)

    use ffdev_topology
    use ffdev_geometry
    use ffdev_geometry_utils

    implicit none
    type(TOPOLOGY)     :: top
    type(GEOMETRY)     :: geo
    ! --------------------------------------------------------------------------

    call ffdev_geometry_utils_comp_dihedrals(.false.,top,geo%trg_crd,geo%crd)

end subroutine ffdev_err_dihedrals_summary

! ------------------------------------------------------------------------------

end module ffdev_err_dihedrals
