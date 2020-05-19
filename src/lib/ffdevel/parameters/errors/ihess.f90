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
    IHessErrorsWeight           = 1.0
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
    integer             :: i,j,q,nihess,ai,aj
    real(DEVDP)         :: err,seterrihess
    real(DEVDP)         :: d0,dt
    ! --------------------------------------------------------------------------

    ! FIXME
    ! error%ihess = 0.0d0

    ! calculate error
    seterrihess = 0.0
    nihess = 0

    do i=1,nsets
        do q=1,sets(i)%top%nbonds
            if( OnlyFFOptIHess) then
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
                seterrihess = seterrihess + sets(i)%geo(j)%weight * err**2
                nihess = nihess + 1
            end do
        end do
    end do

    ! geometry
    if( nihess .gt. 0 ) then
    ! FIXME
       ! error%ihess = sqrt(seterrihess/real(nihess))
    end if

end subroutine ffdev_err_ihess_error

! ==============================================================================
! subroutine ffdev_err_ihess_summary
! ==============================================================================

subroutine ffdev_err_ihess_summary(top,geo,printsum)

    use ffdev_topology
    use ffdev_geometry
    use ffdev_geometry_utils

    implicit none
    type(TOPOLOGY)  :: top
    type(GEOMETRY)  :: geo
    logical         :: printsum
    ! --------------------------------------------------------------------------

    if( .not. geo%trg_crd_optimized ) return

    if( printsum .eqv. .false. ) then
        printsum = .true.
        return
    end if
    ! FIXME
    ! call ffdev_geometry_utils_comp_ihess(.false.,top,geo%trg_crd,geo%crd)

end subroutine ffdev_err_ihess_summary

! ------------------------------------------------------------------------------

end module ffdev_err_ihess
