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

module ffdev_err_impropers

use ffdev_constants

contains

! ==============================================================================
! subroutine ffdev_err_impropers_init
! ==============================================================================

subroutine ffdev_err_impropers_init

    use ffdev_err_impropers_dat

    implicit none
    ! --------------------------------------------------------------------------

    EnableImpropersError            = .false.
    PrintImpropersErrorSummary      = .false.
    ImpropersErrorWeight            = DEV_D2R
    ImpropersErrorLockToPhase       = .true.

end subroutine ffdev_err_impropers_init

! ==============================================================================
! subroutine ffdev_err_impropers_error
! ==============================================================================

subroutine ffdev_err_impropers_error(error)

    use ffdev_targetset_dat
    use ffdev_utils   
    use ffdev_geometry
    use ffdev_errors_dat
    use ffdev_err_impropers_dat

    implicit none
    type(FFERROR_TYPE)  :: error
    ! --------------------------------------------
    integer             :: i,j,q,nimpropers,ai,aj,ak,al,idt
    real(DEVDP)         :: err,seterrimpropers
    real(DEVDP)         :: d0,dt
    ! --------------------------------------------------------------------------

    error%impropers = 0.0

    seterrimpropers = 0.0
    nimpropers = 0

    do i=1,nsets
        do j=1,sets(i)%ngeos

            if( sets(i)%geo(j)%trg_crd_optimized ) then
                do q=1,sets(i)%top%nimpropers
                    ai = sets(i)%top%impropers(q)%ai
                    aj = sets(i)%top%impropers(q)%aj
                    ak = sets(i)%top%impropers(q)%ak
                    al = sets(i)%top%impropers(q)%al
                    d0 = ffdev_geometry_get_improper(sets(i)%geo(j)%crd,ai,aj,ak,al) * DEV_R2D
                    if( ImpropersErrorLockToPhase ) then
                        idt = sets(i)%top%impropers(q)%dt
                        dt = sets(i)%top%improper_types(idt)%g
                    else
                        dt = ffdev_geometry_get_improper(sets(i)%geo(j)%trg_crd,ai,aj,ak,al) * DEV_R2D
                    end if
                    nimpropers = nimpropers + 1
                    err = ffdev_geometry_get_dihedral_deviation(d0,dt)
                    seterrimpropers = seterrimpropers + sets(i)%geo(j)%weight * err**2
                end do
            end if
        end do
    end do

    if( nimpropers .gt. 0 ) then
        error%impropers = sqrt(seterrimpropers/real(nimpropers))
    end if 

end subroutine ffdev_err_impropers_error

! ==============================================================================
! subroutine ffdev_err_impropers_summary
! ==============================================================================

subroutine ffdev_err_impropers_summary(top,geo)

    use ffdev_topology
    use ffdev_geometry
    use ffdev_geometry_utils
    use ffdev_err_impropers_dat

    implicit none
    type(TOPOLOGY)     :: top
    type(GEOMETRY)     :: geo
    ! --------------------------------------------------------------------------

    call ffdev_geometry_utils_comp_impropers(.false.,top,geo%trg_crd,geo%crd,ImpropersErrorLockToPhase)

end subroutine ffdev_err_impropers_summary

! ------------------------------------------------------------------------------

end module ffdev_err_impropers