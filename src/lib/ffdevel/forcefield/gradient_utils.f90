! ==============================================================================
! This file is part of FFDevel.
!    Copyright (C) 2013 Petr Kulhanek, kulhanek@chemi.muni.cz
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

module ffdev_gradient_utils

use ffdev_geometry_dat
use ffdev_constants
use ffdev_variables

contains

! ==============================================================================
! subroutine ffdev_gradient_allocate
! ==============================================================================

subroutine ffdev_gradient_allocate(geo)

    use ffdev_geometry
    use ffdev_utils

    implicit none
    type(GEOMETRY)  :: geo
    ! --------------------------------------------
    integer         :: alloc_stat
    ! --------------------------------------------------------------------------

    if( geo%natoms .le. 0 ) then
        call ffdev_utils_exit(DEV_ERR,1,'Number of atoms has to be greater than zero for gradient!')
    end if

    allocate(geo%grd(3,geo%natoms), stat = alloc_stat)
    if( alloc_stat .ne. 0 ) then
        call ffdev_utils_exit(DEV_ERR,1,'Unable to allocate temporary array for gradient!')
    end if

    geo%grd(:,:) = 0.0d0

end subroutine ffdev_gradient_allocate

! ==============================================================================
! function ffdev_gradient_rmsg
! ==============================================================================

real(DEVDP) function ffdev_gradient_rmsg(geo,maxgrad,maxatom)

    use ffdev_geometry

    implicit none
    type(GEOMETRY)  :: geo
    real(DEVDP)     :: maxgrad
    integer         :: maxatom
    ! ------------------------------
    integer         :: i
    real(DEVDP)     :: norm
    !------------------------------------------------------------------------------

    maxgrad = 0.0
    ffdev_gradient_rmsg = 0.0
    maxatom = 0

    if( geo%natoms .eq. 0 ) return

    do i=1,geo%natoms
        norm = geo%grd(1,i)**2 + geo%grd(2,i)**2 + geo%grd(3,i)**2
        if( abs(geo%grd(1,i)) > abs(maxgrad) ) then
            maxgrad = geo%grd(1,i)
            maxatom = i
        end if
        if( abs(geo%grd(2,i)) > abs(maxgrad) ) then
            maxgrad = geo%grd(2,i)
            maxatom = i
        end if
        if( abs(geo%grd(3,i)) > abs(maxgrad) ) then
            maxgrad = geo%grd(3,i)
            maxatom = i
        end if
        ffdev_gradient_rmsg = ffdev_gradient_rmsg + norm
    end do

    ffdev_gradient_rmsg = sqrt(ffdev_gradient_rmsg/real(3*geo%natoms))

    return

end function ffdev_gradient_rmsg

! ==============================================================================
! function ffdev_gradient_rmsg
! ==============================================================================

real(DEVDP) function ffdev_gradient_rmsg_only(geo)

    use ffdev_geometry

    implicit none
    type(GEOMETRY)  :: geo
    ! ------------------------------
    integer         :: i
    real(DEVDP)     :: norm
    !------------------------------------------------------------------------------

    ffdev_gradient_rmsg_only = 0.0

    if( geo%natoms .eq. 0 ) return

    do i=1,geo%natoms
        norm = geo%grd(1,i)**2 + geo%grd(2,i)**2 + geo%grd(3,i)**2
        ffdev_gradient_rmsg_only = ffdev_gradient_rmsg_only + norm
    end do

    ffdev_gradient_rmsg_only = sqrt(ffdev_gradient_rmsg_only/real(3*geo%natoms))

    return

end function ffdev_gradient_rmsg_only

! ==============================================================================
! subroutine ffdev_gradient_print
! ==============================================================================

subroutine ffdev_gradient_print(fout,top,geo)

    use ffdev_topology
    use ffdev_geometry

    implicit none
    integer         :: fout
    type(TOPOLOGY)  :: top
    type(GEOMETRY)  :: geo
    ! --------------------------------------------
    integer         :: i
    real(DEVDP)     :: stotalnorm
    ! --------------------------------------------------------------------------

    write(DEV_OUT,5) geo%natoms
    write(DEV_OUT,10)
    write(DEV_OUT,20)

    stotalnorm = 0.0d0

    do i=1,geo%natoms
        stotalnorm = stotalnorm + geo%grd(1,i)**2 + geo%grd(2,i)**2 + geo%grd(3,i)**2
        write(fout,30) i, adjustl(top%atoms(i)%name), adjustl(top%atom_types(top%atoms(i)%typeid)%name), &
                          top%atoms(i)%residx, adjustl(top%atoms(i)%resname), &
                          geo%grd(1,i), geo%grd(2,i), geo%grd(3,i)
    end do

    if( geo%natoms .gt. 0 ) then
        stotalnorm = sqrt(stotalnorm/real(3*geo%natoms))
    end if

    write(fout,40)
    write(fout,50) stotalnorm

 5 format('# FFDEVEL GRD ',I8)
10 format('# Idx  Name Type ResI ResN         Gx                Gy                Gz       ')
20 format('# ---- ---- ---- ---- ---- ----------------- ----------------- -----------------')
30 format(I6,1X,A4,1X,A4,1X,I4,1X,A4,1X,F17.8,1X,F17.8,1X,F17.8)
40 format('# ------------------------------------------------------------------------------')
50 format('# Root mean square deviation from zero =                       ',F17.8)
end subroutine ffdev_gradient_print

! ==============================================================================
! subroutine ffdev_gradient_print_notop
! ==============================================================================

subroutine ffdev_gradient_print_notop(fout,geo)

    use smf_periodic_table_dat
    use ffdev_geometry

    implicit none
    integer         :: fout
    type(GEOMETRY)  :: geo
    ! --------------------------------------------
    integer         :: i
    real(DEVDP)     :: stotalnorm
    ! --------------------------------------------------------------------------

    write(DEV_OUT,5) geo%natoms
    write(DEV_OUT,10)
    write(DEV_OUT,20)

    stotalnorm = 0.0d0

    do i=1,geo%natoms
        stotalnorm = stotalnorm + geo%grd(1,i)**2 + geo%grd(2,i)**2 + geo%grd(3,i)**2
        write(fout,30) i, geo%z(i), pt_symbols(geo%z(i)), &
                          geo%grd(1,i), geo%grd(2,i), geo%grd(3,i)
    end do

    if( geo%natoms .gt. 0 ) then
        stotalnorm = sqrt(stotalnorm/real(3*geo%natoms))
    end if

    write(fout,40)
    write(fout,50) stotalnorm

 5 format('# FFDEVEL GRD ',I8)
10 format('# Idx  Z  Sym         Gx                Gy                Gz       ')
20 format('# ---- -- --- ----------------- ----------------- -----------------')
30 format(I6,1X,I2,1X,A3,1X,F17.8,1X,F17.8,1X,F17.8)
40 format('# -----------------------------------------------------------------')
50 format('# Root mean square deviation from zero =          ',F17.8)

end subroutine ffdev_gradient_print_notop

!===============================================================================
!-------------------------------------------------------------------------------
!===============================================================================

logical function ffdev_gradient_test(geo1,geo2,alarm_treshold)

    implicit none
    type(GEOMETRY)  :: geo1
    type(GEOMETRY)  :: geo2
    real(DEVDP)     :: alarm_treshold
    ! ------------------------------
    integer         :: i
    real(DEVDP)     :: diff(3)
    !---------------------------------------------------------------------------

    ffdev_gradient_test = .true.

    if( geo1%natoms .ne. geo2%natoms ) then
        ffdev_gradient_test = .false.
        return
    end if

    do i=1,geo1%natoms
        diff(:) = geo1%grd(:,i) - geo2%grd(:,i)
        if( abs(diff(1)) .gt. alarm_treshold .or. &
            abs(diff(2)) .gt. alarm_treshold .or. &
            abs(diff(3)) .gt. alarm_treshold ) then
            ffdev_gradient_test = .false.
            return
        end if
    end do

end function ffdev_gradient_test

! ------------------------------------------------------------------------------

end module ffdev_gradient_utils
